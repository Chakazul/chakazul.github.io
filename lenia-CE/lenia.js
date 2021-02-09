
var worldSize = 400;
var controlboxWidth = 400;
var controlboxHeight = 400;
var controlGridWidth = 24;
var controlGridHeight = 24;

var N = 64;
var CN = 1;
var KN = 1;
var pixelSize = Math.ceil(worldSize / N);
var updatesPerFrame = 1;
var updateDelay = 0.5;
var worldZoom = 1;
var displayArray = null;
var displayMin = 0;
var displayMax = 1;
var displayN = "CN";
var pattern = null;
var patternID = 0;

var prevColor = 0;
var prevDisplay = 0;
var colormap = d3.interpolateRdYlGn;
var colormapRGB;

var gen = 0;
var time = 0.0;
var timer = null;

// world

var world = d3.selectAll("#cxpbox_lenia_display").append("canvas")
	.attr("width",worldSize)
	.attr("height",worldSize)
    .attr("class","explorable_display");

var canvas = world.node().getContext("2d");

// controls

var controls = d3.selectAll("#cxpbox_lenia_controls").append("svg")
    .attr("width",controlboxWidth)
    .attr("height",controlboxHeight)
    .attr("class","explorable_widgets");	

var g = widget.grid(controlboxWidth, controlboxHeight, controlGridWidth, controlGridHeight);

var playBlock = g.block({x0:5,y0:21,width:0,height:0});
var buttonBlock = g.block({x0:3,y0:16.5,width:4,height:0}).Nx(2);
var radioBlock = g.block({x0:2.5,y0:-0.5,width:0,height:14});
var sliderBlock = g.block({x0:11,y0:14,width:12,height:8}).Ny(4);
var plotBlock = g.block({x0:11,y0:11,width:12,height:1});
var plotBlock2 = g.block({x0:11,y0:11.6,width:12,height:1});
var plotBlock3 = g.block({x0:11,y0:12.2,width:12,height:1});
var toggleBlock = g.block({x0:11,y0:8,width:8,height:0}).Nx(2);
var radioBlock2 = g.block({x0:12,y0:-0.5,width:0,height:6});

//actions: play back pause reload record capture rewind stop
var playButton = { id:"b1", name:"", actions: ["play","stop"], value: 0};
var backButton = { id:"b2", name:"", actions: ["reload"], value: 0};
var stepButton = { id:"b3", name:"step", actions: ["pause"], value: 0};
var playWidget = [
	widget.button(playButton).size(g.x(5)).symbolSize(0.6*g.x(5)).update(runPause)
];
var buttonWidgets = [
	widget.button(backButton).update(() => initPattern(-1)).size(g.x(3)).symbolSize(0.6*g.x(3)),
	widget.button(stepButton).update(oneStep).size(g.x(3)).symbolSize(0.6*g.x(3))
];

var patternRadio = {id: "patterns", name: "patterns", choices: pattern_list.map(d => d.name), value:0};
var displayRadio = {id: "display", name: "display", choices: ["world", "potential", "growth", "kernel"], value:0};
var radioWidgets = [
    widget.radio(patternRadio).size(radioBlock.h()).update(selectPattern)
];
var radioWidgets2 = [
    widget.radio(displayRadio).size(radioBlock2.h()).update(switchDisplay).shape("round")
];

var toggles = [
	{id:"t1", name: "Orli's switch",  value: false},
	{id:"t2", name: "hi-res",  value: false}
];
var toggleWidgets = [
	widget.toggle(toggles[0]).update(switchColor).label("right").size(10),
	widget.toggle(toggles[1]).update(switchRes).label("right").size(10)
];

var slider_m = {id:"m", name: "mu", range: [0.01, 0.7], value: 0};
var slider_s = {id:"s", name: "sigma", range: [0.001, 0.07], value: 0};
var slider_Z = {id:"Z", name: "space scale", range: [20, 1], value: 10};
var slider_T = {id:"T", name: "time scale", range: [20, 1], value: 10};
ww = sliderBlock.w();
vv = true;
var sliderWidget_m = widget.slider(slider_m).width(ww).trackSize(8).handleSize(10).update(setParams).showValue(vv);
var sliderWidget_s = widget.slider(slider_s).width(ww).trackSize(8).handleSize(10).update(setParams).showValue(vv);
var sliderWidget_Z = widget.slider(slider_Z).width(ww).trackSize(8).handleSize(10).update(setZoom).showValue(vv);
var sliderWidget_T = widget.slider(slider_T).width(ww).trackSize(8).handleSize(10).update(setParams).showValue(vv);
var sliderWidgets1 = [sliderWidget_s, sliderWidget_m];
var sliderWidgets2 = [sliderWidget_T, sliderWidget_Z];

var pb = controls.selectAll(".button .play").data(playWidget).enter().append(widget.buttonElement)
	.attr("transform",(d,i) => "translate("+playBlock.x(0)+","+playBlock.y(0)+")");	

var bu = controls.selectAll(".button .others").data(buttonWidgets).enter().append(widget.buttonElement)
	.attr("transform",(d,i) => "translate("+buttonBlock.x(i)+","+buttonBlock.y(0)+")");	

var rad = controls.selectAll(".radio").data(radioWidgets).enter().append(widget.radioElement)
	.attr("transform",(d,i) => "translate("+radioBlock.x(0)+","+radioBlock.y(0)+")");	

var rad2 = controls.selectAll(".radio .two").data(radioWidgets2).enter().append(widget.radioElement)
	.attr("transform",(d,i) => "translate("+radioBlock2.x(0)+","+radioBlock2.y(0)+")");	

var tg = controls.selectAll(".toggle").data(toggleWidgets).enter().append(widget.toggleElement)
	.attr("transform",(d,i) => "translate("+toggleBlock.x(i)+","+toggleBlock.y(0)+")");	

var sl = controls.selectAll(".slider").data(sliderWidgets1).enter().append(widget.sliderElement)
	.attr("transform",(d,i) => "translate("+sliderBlock.x(0)+","+sliderBlock.y(i)+")");

var sl2 = controls.selectAll(".slider .two").data(sliderWidgets2).enter().append(widget.sliderElement)
	.attr("transform",(d,i) => "translate("+sliderBlock.x(0)+","+sliderBlock.y(i+2)+")");

var colors1 = controls.append("g").attr("class","bars")
	.attr("transform",(d,i) => "translate("+plotBlock.x(0)+","+plotBlock.y(0)+")");

var colors2 = controls.append("g").attr("class","bars")
	.attr("transform",(d,i) => "translate("+plotBlock2.x(0)+","+plotBlock2.y(0)+")");

var colors3 = controls.append("g").attr("class","bars")
	.attr("transform",(d,i) => "translate("+plotBlock3.x(0)+","+plotBlock3.y(0)+")");

var bar1 = colors1.selectAll(".bars")
    .data(d3.range(plotBlock.w()), d => d)
    .enter().append("rect")
    .attr("class", "bars")
    .attr("x", (d,i) => i)
    .attr("y", 0)
    .attr("height", 10)
    .attr("width", 1)
    .style("fill", (d,i) => colormap(d/plotBlock.w()));

var bar2 = colors2.selectAll(".bars")
    .data(d3.range(plotBlock.w()), d => d)
    .enter().append("rect")
    .attr("class", "bars")
    .attr("x", (d,i) => i)
    .attr("y", 0)
    .attr("height", 0)
    .attr("width", 1);

var bar3 = colors3.selectAll(".bars")
    .data(d3.range(plotBlock.w()), d => d)
    .enter().append("rect")
    .attr("class", "bars")
    .attr("x", (d,i) => i)
    .attr("y", 0)
    .attr("height", 0)
    .attr("width", 1);

var lh = colors1.selectAll(".plot").data(["low","high"]).enter().append("text")
	.text(d => d)
	.attr("transform", (d,i) => "translate("+(i*plotBlock.w())+",25)")
	.attr("class","plot")
	.style("text-anchor","middle");

// maths

const PRECISION = 1000000;
function round(x) { return Math.round(x * PRECISION) / PRECISION; }
function mod(x, n) { return ((x % n) + n) % n; }

function fft1d(dir, re1, im1) {
	/* Do the bit reversal */
	var n = re1.length, m = round(Math.log2(n)), n2 = n >> 1, j1 = 0;
	for (var j=0; j<n-1; j++) {
		if (j < j1) {
			var tmp = re1[j]; re1[j] = re1[j1]; re1[j1] = tmp;
			tmp = im1[j]; im1[j] = im1[j1]; im1[j1] = tmp;
		}
		var j2 = n2;
		while (j2 <= j1) {
			j1 -= j2;
			j2 >>= 1;
		}
		j1 += j2;
	}
	
	/* Compute the FFT */
	var c1 = -1.0, c2 = 0.0, l2 = 1;
	for (var l=0; l<m; l++) {
		var l1 = l2;
		l2 <<= 1;
		var u1 = 1.0, u2 = 0.0;
		for (var i=0; i<l1; i++) {
			for (var j=i; j<n; j+=l2) {
				var j2 = j + l1;
				var t1 = u1 * re1[j2] - u2 * im1[j2];
				var t2 = u1 * im1[j2] + u2 * re1[j2];
				re1[j2] = re1[j] - t1;
				im1[j2] = im1[j] - t2;
				re1[j] += t1;
				im1[j] += t2;
			}
			var z = u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = z;
		}
		c2 = Math.sqrt((1.0 - c1) / 2.0);
		if (dir == 1)
			c2 = -c2;
		c1 = Math.sqrt((1.0 + c1) / 2.0);
	}
	
	/* Scaling for forward transform */
	if (dir == -1) {
		var scale_f = 1.0 / n;		
		for (var j=0; j<n; j++) {
			re1[j] *= scale_f;
			im1[j] *= scale_f;
		}
	}
}

function fft2d(dir, re2, im2) {
	var n = re2.length;
	for (var i=0; i<n; i++)
		fft1d(dir, re2[i], im2[i]);
    for (var i=0; i<n; i++) {
        for (var j=0; j<i; j++) {
    		var tmp = re2[i][j]; re2[i][j] = re2[j][i]; re2[j][i] = tmp;
    	}
    }
	for (var i=0; i<n; i++) {
        for (var j=0; j<i; j++) {
    		var tmp = im2[i][j]; im2[i][j] = im2[j][i]; im2[j][i] = tmp;
    	}
    }
	for (var i=0; i<n; i++)
		fft1d(dir, re2[i], im2[i]);
}

function matrixMult(ar, ai, br, bi, cr, ci) {
	var n = ar.length;
	for (var i=0; i<n; i++) {
		var ar_i = ar[i], ai_i = ai[i];
		var br_i = br[i], bi_i = bi[i];
		var cr_i = cr[i], ci_i = ci[i];
		for (var j=0; j<n; j++) {
			var a = ar_i[j]; var b = ai_i[j];
			var c = br_i[j]; var d = bi_i[j];
			var t = a * (c + d);
			cr_i[j] = t - d*(a+b);
			ci_i[j] = t + c*(b-a);
		}
	}
}

// core functions

function deltaFunc(n, m, s) {
	var r = n - m;
    return Math.exp(-r*r/ (2*s*s) ) * 2 - 1;
}

function kernelCoreFunc(r) {
	return Math.exp(4 - 1/r/(1-r));
}

function kernelFunc(r, b) {
	if (r>=1) return 0;
	var R = r * b.length;
	return kernelCoreFunc(R % 1) * b[Math.floor(R)];
}

function calcKernel(zoom) {
    for (var k=0; k<KN; k++) {
    	var weight = 0.0;
    	var kernel = pattern.kernels[k];
    	for (var y=0; y<N; y++) {
    		for (var x=0; x<N; x++) {
    			var yy = ((y + N/2) % N) - N/2;
    			var xx = ((x + N/2) % N) - N/2;
    			var r = Math.sqrt(xx*xx + yy*yy) / pattern.defaultR / zoom / worldZoom / kernel.r;
    			var v = kernelFunc(r, kernel.b);
    			weight += v;
    			Kr[k][y][x] = v;
    			yy = N - ((y + N/2) % N) - 1;
    			xx = ((x + N/2) % N);
    			K[k][yy][xx] = v;
    		}
    	}

    	fft2d(1, Kr[k], Ki[k]);
    	for (var y=0; y<N; y++) {
    		for (var x=0; x<N; x++) {
    			Kr[k][y][x] /= weight;
    			Ki[k][y][x] /= weight;
    		}
    	}
    }
}

function update() {
    for (var c=0; c<CN; c++)
    	for (var y=0; y<N; y++)
    		for (var x=0; x<N; x++)
    			Ao[c][y][x] = Ar[c][y][x];
    for (var c=0; c<CN; c++) {
        for (var y=0; y<N; y++) {
    		Ai[c][y].fill(0);
    		G[c][y].fill(0);
        }
    }

	// f * g = F-1( F(f) dot F(g) )
    for (var c=0; c<CN; c++)
        fft2d(1, Ar[c], Ai[c]);
    for (var k=0; k<KN; k++) {
        var kernel = pattern.kernels[k];
        matrixMult(Ar[kernel.c0], Ai[kernel.c0], Kr[k], Ki[k], Ur[k], Ui[k]);
    	fft2d(-1, Ur[k], Ui[k]);
    }

    for (var k=0; k<KN; k++) {
        var kernel = pattern.kernels[k];
        var m = KN == 1 ? global_m : kernel.m;
        var s = KN == 1 ? global_s : kernel.s;
    	for (var y=0; y<N; y++) {
    		for (var x=0; x<N; x++) {
    			var u = Ur[k][y][x];
    			G[kernel.c1][y][x] += deltaFunc(u, m, s) * kernel.h;
    		}
    	}
    }

    for (var c=0; c<CN; c++) {
    	for (var y=0; y<N; y++) {
    		for (var x=0; x<N; x++) {
                var g = G[c][y][x];
        		var v = Ao[c][y][x] + g / global_T * updateDelay;
        		if (v<0) v = 0; else if (v>1) v = 1;
        		Ar[c][y][x] = v;
            }
        }
    }

	gen++;
	time = round(time + round(1/global_T));
}

// functions

function interpolateTurbo(t) {
    t = Math.max(0, Math.min(1, t));
    return "rgb("
        + Math.max(0, Math.min(255, Math.round(34.61 + t * (1172.33 - t * (10793.56 - t * (33300.12 - t * (38394.49 - t * 14825.05))))))) + ", "
        + Math.max(0, Math.min(255, Math.round(23.31 + t * (557.33 + t * (1225.33 - t * (3574.96 - t * (1073.77 + t * 707.56))))))) + ", "
        + Math.max(0, Math.min(255, Math.round(27.2 + t * (3211.1 - t * (15327.97 - t * (27814 - t * (22569.18 - t * 6838.66)))))))
        + ")";
}

function initWorld(N0, frame0, delay0, zoom0) {
	N = N0;
	pixelSize = Math.ceil(worldSize / N);
	updatesPerFrame = frame0;
	updateDelay = delay0;
	worldZoom = zoom0;

	Ar = initArray(3);
	Ai = initArray(3);
	Ao = initArray(3);
	K  = initArray(16);
	Kr = initArray(16);
	Ki = initArray(16);
	Ur = initArray(16);
	Ui = initArray(16);
	G  = initArray(3);
	displayArray = Ar;
	displayMin = 0;
	displayMax = 1;
    displayN = "CN";

	xScale = d3.scaleLinear().domain([0,N]).range([0,worldSize]);
	yScale = d3.scaleLinear().domain([0,N]).range([0,worldSize]);
}

function initArray(num) {
    var arr = [];
	for (var j=0; j<num; j++) {
        var ch = [];
        for (var i=0; i<N; i++)
            ch.push(new Array(N).fill(0));
        arr.push(ch);
    }
	return arr;
}

function switchRes(d) {
	if (d.value == 0) {
		initWorld(64, 1, 0.5, 1);
	} else {
		//initWorld(256, 1, 1, 4);
		initWorld(128, 1, 1, 2);
	}
    switchDisplay();  // need to re-attach displayArray after re-init arrays
	initPattern();
}

function switchColor(d) {
    var i = d ? d.value : prevColor;
	if (i == 0) {
		colormap = d3.interpolateRdYlGn;
        var gamma = 0.8;
        colormapRGB = [
            d3.interpolateRgb.gamma(gamma)("black","rgb(100%,0%,0%)"), 
            d3.interpolateRgb.gamma(gamma)("black","rgb(0%,100%,0%)"), 
            d3.interpolateRgb.gamma(gamma)("black","rgb(0%,0%,100%)")
        ];
	} else {
		//colormap = d3.interpolateRainbow;
        colormap = interpolateTurbo;
        var gamma = 0.7;
        colormapRGB = [
            d3.interpolateRgb.gamma(gamma)("black","rgb(75%,38%,0%)"), 
            d3.interpolateRgb.gamma(gamma)("black","rgb(0%,75%,38%)"), 
            d3.interpolateRgb.gamma(gamma)("black","rgb(38%,0%,75%)")
        ]; // OTV = Orange/Turquoise/Violet
	}

    if (CN == 1) {
    	colors1.selectAll(".bars").attr("height", 10).style("fill", (d,i) => colormap(d/plotBlock.w()));
        colors2.selectAll(".bars").attr("height", 0);
        colors3.selectAll(".bars").attr("height", 0);
    } else {
    	colors1.selectAll(".bars").attr("height", 10).style("fill", (d,i) => colormapRGB[2](d/plotBlock.w()));
    	colors2.selectAll(".bars").attr("height", 10).style("fill", (d,i) => colormapRGB[1](d/plotBlock.w()));
    	colors3.selectAll(".bars").attr("height", 10).style("fill", (d,i) => colormapRGB[0](d/plotBlock.w()));
    }
	draw();

    prevColor = i;
}

function switchDisplay(d) {
    var i = d ? d.value() : prevDisplay;
	switch (i) {
		case 0: displayArray = Ar; displayN = "CN"; displayMin = 0;  displayMax = 1; break;
		case 1: displayArray = Ur; displayN = "KN"; displayMin = 0;  displayMax = global_m*2; break;
		case 2: displayArray = G;  displayN = "CN"; displayMin = -1; displayMax = 1; break;
		case 3: displayArray = K;  displayN = "KN"; displayMin = 0;  displayMax = 1; break;
	}
	draw();
    prevDisplay = i;
}

function setParams(d) {
    if (KN == 1) {
    	global_m = slider_m.value;
    	if (displayRadio.value == 1) displayMax = global_m * 2;
    	global_s = slider_s.value;
    }
	global_T = slider_T.value;
}
function setZoom(d) {
	global_Z = slider_Z.value;
    calcKernel(global_Z/10);
}

function run() {
    for (var i=0; i<updatesPerFrame; i++)
        update();
    draw();
}

function runPause(d) {
    d.value == 1 ? timer = d3.timer(run, 1) : timer.stop();
}

function oneStep() {
	if (playButton.value == 1)
		playWidget[0].click();
	run();
}

function selectPattern(d) {
	initPattern(d.value());
}

function draw() {
	//canvas.fillStyle = 0;
	//canvas.fillRect(0, 0, worldSize, worldSize);
    if (KN == 1) {
    	for (var y=0; y<N; y++) {
            for (var x=0; x<N; x++) {
        		var v = displayArray[0][y][x];
        		v = (v - displayMin) / (displayMax - displayMin);
                canvas.fillStyle = colormap(v);
                canvas.fillRect(xScale(x), yScale(y), pixelSize, pixelSize);
            }
        }
    } else if (displayN == "CN") {
    	for (var y=0; y<N; y++) {
            for (var x=0; x<N; x++) {
                var color = d3.rgb(0, 255*0.2, 0);
            	for (var c=0; c<CN; c++) {
            		var v = displayArray[c][y][x];
            		v = (v - displayMin) / (displayMax - displayMin);
                    color0 = d3.rgb(colormapRGB[c](v));
                    color.r += color0.r;
                    color.g += color0.g;
                    color.b += color0.b;
                }
                canvas.fillStyle = color;
                canvas.fillRect(xScale(x), yScale(y), pixelSize, pixelSize);
            }
        }
    } else if (displayN == "KN") {
    	for (var y=0; y<N; y++) {
            for (var x=0; x<N; x++) {
                var color = d3.rgb(0, 255*0.2, 0);
            	for (var k=0; k<KN; k++) {
                    var kernel = pattern.kernels[k];
            		var v = displayArray[k][y][x] / 2;
                    color0 = d3.rgb(colormapRGB[kernel.c0](v));
                    color.r += color0.r;
                    color.g += color0.g;
                    color.b += color0.b;
                }
                canvas.fillStyle = color;
                canvas.fillRect(xScale(x), yScale(y), pixelSize, pixelSize);
            }
        }
    }
}

function clearWorld() {
	for (var c=0; c<CN; c++)
    	for (var y=0; y<N; y++)
            Ar[c][y].fill(0);
}

function placePattern(p, shifty, shiftx, angle, zoom) {
    var scale = p.defaultR / p.R * zoom * worldZoom;
	var sin = Math.sin(angle / 180 * Math.PI);
	var cos = Math.cos(angle / 180 * Math.PI);
    var ch = p.cells[0].length;
    var cw = p.cells[0][0].length;
	var h = (Math.abs(ch*cos) + Math.abs(cw*sin) + 1) * scale - 1;
	var w = (Math.abs(cw*cos) + Math.abs(ch*sin) + 1) * scale - 1;
    var y0 = Math.floor(N/2 - h/2 + shifty);
    var x0 = Math.floor(N/2 - w/2 + shiftx);
	for (var c=0; c<CN; c++) {
    	for (var y=0; y<h; y++) {
            for (var x=0; x<w; x++) {
        		var cy = Math.round( (- (x-w/2)*sin + (y-h/2)*cos) / scale + ch/2 );
        		var cx = Math.round( (+ (x-w/2)*cos + (y-h/2)*sin) / scale + cw/2 );
        		var v = (cy>=0 && cx>=0 && cy<ch && cx<cw) ? p.cells[c][cy][cx] : 0;
        		Ar[c][(y0+y+N)%N][(x0+x+N)%N] = v;
        	}
        }
    }
}

function randomPlacePattern(p) {
	var chance = Math.random();
    var zoom = 1;
	if (pattern.name.includes("(s)") && chance>0.7) {
		zoom = 0.7;
		placePattern(p, -N/4, 0, Math.random()*360, zoom);
		placePattern(p, +N/4, +N/4, Math.random()*360, zoom);
		placePattern(p, +N/4, -N/4, Math.random()*360, zoom);
	} else if (pattern.name.includes("(s)") && chance>0.3) {
		zoom = 0.8;
		placePattern(p, +N/4, +N/4, Math.random()*360, zoom);
		placePattern(p, -N/4, -N/4, Math.random()*360, zoom);
	} else {
		placePattern(p, 0, 0, Math.random()*360, zoom);
	}

	calcKernel(zoom);
    sliderWidget_Z.click(zoom*10);
}

function showSliders(isShow){
	if (isShow){
		sl.transition().style("opacity", (d,i) => 1)
            .select(".track-overlay")
            .style("pointer-events","all");
	} else {
		sl.transition().style("opacity", (d,i) => i>3 ? 1 : 0)
            .select(".track-overlay")
            .style("pointer-events", (d,i) => i>3 ? "all" : "none");
	}
}

function initPattern(id, isInit=false) {
	if (id>=0) {
		patternID = id;
		pattern = pattern_list[id];
        CN = pattern.cells.length;
        KN = pattern.kernels.length;
        switchColor();
	}
    if (KN == 1) {
	    var kernel = pattern.kernels[0];
    	sliderWidget_m.click(kernel.m);
    	sliderWidget_s.click(kernel.s);
        showSliders(true);
    } else {
        showSliders(false);
    }
	sliderWidget_T.click(pattern.T);

    clearWorld();
	if (isInit) {
        var zoom = 1;
        placePattern(pattern, 0, 0, 0, zoom);
    	calcKernel(zoom);
        sliderWidget_Z.click(zoom*10);
    } else {
        randomPlacePattern(pattern);
    }
    
	update();
    draw();
    gen = 0;
    time = 0.0;
}

initWorld(64, 1, 0.5, 1);
initPattern(0, isInit=true);

