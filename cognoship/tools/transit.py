"""Experiment: how early can a ??? stretch's eventual outcome be told from its trajectory?
For each sweep trial (G, R, B) of a morph whose name box goes ??? after the action: features over the
last 20 steps at tau steps into the ???, the eventual named outcome (the trial's label). Predictors,
judged on the trials still ??? at tau:
  prior   the action's majority outcome over its other starts (no watching: our model)
  generic kNN over trials of OTHER actions (a transit detector that knows no action)
  action  kNN over the same action's other starts (an action-specific transit)
  all     kNN over every other trial
Usage: tools/transit.py CODE MORPH (from Cognoship/; results in output/transit/)"""
import sys, json, os, collections
import numpy as np, torch
code, morph = sys.argv[1], sys.argv[2]
D = f'output/models/{code}'
TAUS, W, K = (5, 10, 20, 40, 80), 20, 5
FEATS = ('mass', 'mR', 'mG', 'mB', 'speed', 'omega')
X, Y, ACT, T0, END = [], [], [], [], []
for tag, ch in (('', 'G'), ('-r', 'R'), ('-b', 'B')):
    jp, zp = f'{D}/steer-{morph}{tag}.json', f'{D}/steer-{morph}{tag}-series.npz'
    if not os.path.exists(zp): continue
    runs = json.load(open(jp))['runs']; z = np.load(zp)
    named = z['named']; F = np.stack([np.asarray(z[k], np.float32) for k in FEATS], 2)
    F[:, :, 5] = np.abs(F[:, :, 5])   # (turning either way alike)
    for i, r in enumerate(runs):
        if r['action'] == 'none': continue
        unk = np.where(named[i] == '')[0]
        unk = unk[unk >= 5]   # (after the start)
        if not len(unk): continue
        t0 = int(unk[0]); te = t0
        while te < named.shape[1] and named[i, te] == '': te += 1   # (the stretch's end: named again)
        X.append(F[i]); Y.append(r['label']); ACT.append((ch, r['spot'], r['dur'])); T0.append(t0); END.append(te)
X = np.stack(X); Y = np.array(Y); T0 = np.array(T0); END = np.array(END)
n, L = len(Y), X.shape[1]
mu, sd = X.reshape(-1, len(FEATS)).mean(0), X.reshape(-1, len(FEATS)).std(0) + 1e-9
X = (X - mu) / sd
akey = np.array([hash(a) for a in ACT])
print(f'{code} {morph}: {n} trials with a ??? stretch; its length median {int(np.median(END - T0))} steps (quartiles {np.percentile(END - T0, [25, 75]).astype(int).tolist()})')
print('outcomes:', dict(collections.Counter(Y).most_common(8)))
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
res = {}
for tau in TAUS:
    t = T0 + tau
    live = (t < L) & (END > t)   # (still ??? at tau)
    idx = np.where(live)[0]
    if len(idx) < 30: continue
    win = np.stack([X[i, t[i] - W + 1:t[i] + 1:4].ravel() if t[i] - W + 1 >= 0 else np.zeros(((W + 3) // 4) * len(FEATS)) for i in idx])
    V = torch.tensor(win, device=dev); d = torch.cdist(V, V); d.fill_diagonal_(float('inf'))
    a, yy = akey[idx], Y[idx]
    same = torch.tensor(a[:, None] == a[None, :], device=dev)
    def knn(mask):
        dd = torch.where(mask, d, torch.tensor(float('inf'), device=dev))
        v, j = dd.topk(K, largest=False)
        out = []
        for r in range(len(idx)):
            jj = j[r][torch.isfinite(v[r])].cpu().numpy()
            out.append(collections.Counter(yy[jj]).most_common(1)[0][0] if len(jj) else None)
        return np.array(out, dtype=object)
    pred = {'generic': knn(~same), 'action': knn(same), 'all': knn(torch.ones_like(same))}
    prior = []
    for r in range(len(idx)):
        o = [yy[j] for j in np.where(a == a[r])[0] if j != r]
        prior.append(collections.Counter(o).most_common(1)[0][0] if o else None)
    pred['prior'] = np.array(prior, dtype=object)
    # the uncertain actions (their own outcome under 90% over their starts): where watching could add
    unc = np.array([collections.Counter(yy[a == a[r]]).most_common(1)[0][1] / (a == a[r]).sum() < 0.9 for r in range(len(idx))])
    row = {'still': int(live.sum()), 'still_share': round(live.mean(), 2)}
    for k, p in pred.items():
        row[k] = round(float(np.mean(p == yy)), 2); row[k + '_unc'] = round(float(np.mean(p[unc] == yy[unc])), 2) if unc.any() else None
    row['n_unc'] = int(unc.sum())
    # spread: each trial's nearest same-outcome neighbour against its nearest other-outcome one
    dn = d.clone(); so = torch.tensor(yy[:, None] == yy[None, :], device=dev)
    ds = torch.where(so, dn, torch.tensor(float('inf'), device=dev)).min(1).values
    do = torch.where(~so, dn, torch.tensor(float('inf'), device=dev)).min(1).values
    fin = torch.isfinite(ds) & torch.isfinite(do)
    row['separation'] = round(float((do[fin] / ds[fin]).median()), 2)
    res[tau] = row
    print(f"tau {tau:3d}: still ??? {row['still']} ({row['still_share']:.0%}) | all: prior {row['prior']:.2f} generic {row['generic']:.2f} action {row['action']:.2f} all {row['all']:.2f}"
          f" | uncertain actions ({row['n_unc']}): prior {row['prior_unc']} generic {row['generic_unc']} action {row['action_unc']} all {row['all_unc']} | other/same NN distance {row['separation']}", flush=True)
json.dump(res, open(f'output/transit/{code}-{morph}.json', 'w'), indent=1)
