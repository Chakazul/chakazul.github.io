"""Split the CARL policy (models/agent_direction.onnx) into a FiLM model and a conv core.

The exported policy is a 4-level U-Net whose 14 conv blocks are each FiLM-conditioned on the
4-float `context`: a 3-layer MLP per block turns the context into gamma/beta, applied as
Mul + Add after the conv. On WebGPU every ONNX node is a separate dispatch, and the export spends
~680 of its 694 nodes on things other than the 15 convs: the 14 MLPs (56 Gemm), the shape
arithmetic that reshapes their outputs and sizes the upsamples, and circular padding built from
4 Slice + 2 Concat per conv. This script produces, with no change to the arithmetic:

  agent_direction_film.onnx         context -> film_g{i}, film_b{i} ([1,C,1,1], i = 0..13).
                                    Run once per context change, not per inference.
  agent_direction_core.onnx         state + film_* -> q. Circular padding kept as Slice/Concat,
                                    spatial dims still dynamic (any crop size). 150 ops.
  agent_direction_core_gather.onnx  Same, but each circular pad is two Gathers with constant
                                    indices [n-1, 0..n-1, 0], which fixes the input at 96x96.
                                    94 ops: the fewest dispatches, but Gather is slow on CPU,
                                    so this one is for WebGPU only.

The upsamples' target size (computed from the skip connection's shape in the export) becomes a
constant scale of 2, which is the same thing at every level of this net. Both cores reproduce the
original's q bit for bit (checked on 400 crops from closed-loop rollouts, see PERF.md).

Usage, from this folder:
    uv run --with onnx --with numpy python split_model.py
"""
import os
import numpy as np
import onnx
from onnx import helper, numpy_helper as nh, TensorProto

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(HERE, '..', 'models')
SRC = os.path.join(MODELS, 'agent_direction.onnx')
NET = 96   # the gather core's fixed crop size (CFG.netSize's default)


def load():
    m = onnx.load(SRC)
    g = m.graph
    prod = {o: n for n in g.node for o in n.output}
    cons = {}
    for n in g.node:
        for i in n.input:
            cons.setdefault(i, []).append(n)
    return m, g, prod, cons


def upstream(prod, names, stop=lambda x: False):
    """Nodes that `names` depend on, stopping at graph inputs, initializers or `stop`."""
    seen, out, stack = set(), [], list(names)
    while stack:
        x = stack.pop()
        n = prod.get(x)
        if n is None or id(n) in seen or stop(x):
            continue
        seen.add(id(n))
        out.append(n)
        stack.extend(n.input)
    return out


def film_sites(g, prod, cons):
    """(Mul, Add, gamma Reshape, beta Reshape, channels) for each FiLM block, in graph order."""
    weights = {i.name: i for i in g.initializer}
    sites = []
    for n in g.node:
        if n.op_type != 'Mul':
            continue
        rg = prod.get(n.input[0])
        if not (rg and rg.op_type == 'Reshape' and prod[rg.input[0]].op_type == 'Gemm'):
            continue
        add = cons[n.output[0]][0]
        rb = prod[add.input[1]]
        assert add.op_type == 'Add' and rb.op_type == 'Reshape'
        C = nh.to_array(weights[prod[rg.input[0]].input[1]]).shape[0]
        sites.append((n, add, rg, rb, int(C)))
    assert len(sites) == 14, len(sites)
    return sites


def film_inputs(sites):
    return [helper.make_tensor_value_info(f'film_{t}{i}', TensorProto.FLOAT, [1, s[4], 1, 1])
            for i, s in enumerate(sites) for t in 'gb']


def make_film(m, g, prod, sites):
    nodes, outs, inits, used = [], [], {}, set()
    for i, (_, _, rg, rb, C) in enumerate(sites):
        for tag, r in (('g', rg), ('b', rb)):
            for n in reversed(upstream(prod, [r.input[0]], stop=lambda x: x == 'context')):
                if id(n) not in used:
                    used.add(id(n))
                    nodes.append(n)
            shp = f'film_shape_{C}'
            inits.setdefault(shp, nh.from_array(np.array([1, C, 1, 1], np.int64), shp))
            nodes.append(helper.make_node('Reshape', [r.input[0], shp], [f'film_{tag}{i}']))
            outs.append(helper.make_tensor_value_info(f'film_{tag}{i}', TensorProto.FLOAT, [1, C, 1, 1]))
    need = {i for n in nodes for i in n.input}
    weights = [i for i in g.initializer if i.name in need] + list(inits.values())
    graph = helper.make_graph(nodes, 'film',
        [helper.make_tensor_value_info('context', TensorProto.FLOAT, [1, 4])], outs, weights)
    return helper.make_model(graph, opset_imports=m.opset_import, ir_version=m.ir_version)


def make_core(gather):
    m, g, prod, cons = load()
    sites = film_sites(g, prod, cons)
    for i, (mul, add, _, _, _) in enumerate(sites):
        mul.input[0] = f'film_g{i}'
        add.input[1] = f'film_b{i}'
    scales = nh.from_array(np.array([1, 1, 2, 2], np.float32), 'up_scales')
    for n in g.node:
        if n.op_type == 'Resize':      # (X, roi, scales, sizes) -> (X, '', scales)
            n.input[1] = ''
            n.input[2] = 'up_scales'
            del n.input[3:]

    nodes, extra_inits = list(g.node), {}
    if gather:
        # Spatial size at each conv, from shape inference on the rewired graph at NET x NET.
        probe = helper.make_model(helper.make_graph(nodes, 'probe',
            [helper.make_tensor_value_info('state', TensorProto.FLOAT, [1, 4, NET, NET])] + film_inputs(sites),
            [helper.make_tensor_value_info('q', TensorProto.FLOAT, None)], list(g.initializer) + [scales]),
            opset_imports=m.opset_import, ir_version=m.ir_version)
        shape = {v.name: [d.dim_value for d in v.type.tensor_type.shape.dim]
                 for v in onnx.shape_inference.infer_shapes(probe).graph.value_info}
        shape['state'] = [1, 4, NET, NET]
        for conv in [n for n in nodes if n.op_type == 'Conv']:
            ks = [a for a in conv.attribute if a.name == 'kernel_shape'][0].ints
            if list(ks) != [3, 3]:
                continue
            # conv <- Concat(axis 3) <- Concat(axis 2) whose middle input is the unpadded tensor
            x = prod[prod[conv.input[0]].input[1]].input[1]
            H, W = shape[x][2], shape[x][3]
            for n in (H, W):
                name = f'wrap_idx_{n}'
                # int32 rather than int64: WGSL has no 64-bit integers, so int64 indices cost a
                # conversion on WebGPU. ONNX Gather accepts either.
                extra_inits.setdefault(name, nh.from_array(np.array([n - 1, *range(n), 0], np.int32), name))
            gy = helper.make_node('Gather', [x, f'wrap_idx_{H}'], [x + '_wrapH'], axis=2)
            gx = helper.make_node('Gather', [x + '_wrapH', f'wrap_idx_{W}'], [x + '_wrapHW'], axis=3)
            conv.input[0] = gx.output[0]
            at = next(i for i, n in enumerate(nodes) if n is conv)
            nodes[at:at] = [gy, gx]

    # Keep only what q depends on: drops the MLPs, the shape arithmetic and (gather) the Slices.
    prod = {o: n for n in nodes for o in n.output}
    keep = {id(n) for n in upstream(prod, ['q'])}
    nodes = [n for n in nodes if id(n) in keep]
    need = {i for n in nodes for i in n.input}
    inits = [i for i in g.initializer if i.name in need] + [scales] + list(extra_inits.values())
    hw = [NET, NET] if gather else ['H', 'W']
    graph = helper.make_graph(nodes, 'core',
        [helper.make_tensor_value_info('state', TensorProto.FLOAT, [1, 4, *hw])] + film_inputs(sites),
        [helper.make_tensor_value_info('q', TensorProto.FLOAT, [1, 3, *hw])], inits)
    return helper.make_model(graph, opset_imports=m.opset_import, ir_version=m.ir_version)


def main():
    m, g, prod, cons = load()
    out = {
        'agent_direction_film.onnx': make_film(m, g, prod, film_sites(g, prod, cons)),
        'agent_direction_core.onnx': make_core(gather=False),
        'agent_direction_core_gather.onnx': make_core(gather=True),
    }
    for name, model in out.items():
        onnx.checker.check_model(model)
        onnx.save(model, os.path.join(MODELS, name))
        ops = [n.op_type for n in model.graph.node if n.op_type != 'Constant']
        print(f'{name}: {len(ops)} ops')


if __name__ == '__main__':
    main()
