(function(){

// variables

var worldSize = 400;
var controlboxWidth = 400;
var controlboxHeight = 400;
var controlGridWidth = 24;
var controlGridHeight = 24;

var N = 64;
var CN = 1;
var KN = 1;
var pixelSize = Math.ceil(worldSize / N);
var worldDelay = 0.5;
var worldZoom = 1;
var displayArray = null;
var displayMin = 0;
var displayMax = 1;
var displayN = "CN";
var pattern = null;
var isDiscrete = false;

var global_m = 0.15;
var global_s = 0.015;
var global_R = 1;
var global_T = 10;
var global_P = 50;

var colormap = d3.interpolateViridis;
var colormapRGB;

var gen = 0;
var time = 0.0;
var timer = null;

// world

var display = d3.selectAll("#cxpbox_lenia_display").append("canvas")
	.attr("width",worldSize)
	.attr("height",worldSize)
    .attr("class","explorable_display");

var canvas = display.node().getContext("2d");

// controls

var controls = d3.selectAll("#cxpbox_lenia_controls").append("svg")
    .attr("width",controlboxWidth)
    .attr("height",controlboxHeight)
    .attr("class","explorable_widgets");	

var g = widget.grid(controlboxWidth, controlboxHeight, controlGridWidth, controlGridHeight);

// control position

var playBlock = g.block({x0:5, y0:21, width:0, height:0});
var buttonBlock = g.block({x0:3, y0:16.5, width:4, height:0}).Nx(2);
var radioBlock = g.block({x0:2.5, y0:-0.5, width:0, height:14});
var sliderBlock = g.block({x0:11, y0:12, width:12, height:10}).Ny(5);
var plotBlock = g.block({x0:11, y0:10, width:12, height:1});
var plotBlock2 = g.block({x0:11, y0:10.6, width:12, height:1});
var plotBlock3 = g.block({x0:11, y0:11.2, width:12, height:1});
var radioBlock2 = g.block({x0:12, y0:-0.5, width:0, height:7});
var toggleBlock = g.block({x0:20, y0:3, width:4, height:3}).Ny(2);
//var toggleBlock = g.block({x0:17.5, y0:6.5, width:4, height:3}).Nx(2);
//var legendBlock = g.block({x0:16, y0:1, width:4, height:2}).Ny(3);

// control data & widget

//button actions: play back pause reload record capture rewind stop
var playButton = { id:"b1", name:"", actions: ["play","stop"], value: 0};
//var reloadButton = { id:"b2", name:"", actions: ["reload"], value: 0};
var resetButton = { id:"b3", name:"", actions: ["back"], value: 0};
var stepButton = { id:"b4", name:"step", actions: ["pause"], value: 0};
var playWidget = [
	widget.button(playButton).size(g.x(5)).symbolSize(0.6*g.x(5)).update(runPause)
];
var buttonWidgets = [
	//widget.button(reloadButton).update(initPattern).size(g.x(3)).symbolSize(0.6*g.x(3)),
	widget.button(resetButton).update(initCells).size(g.x(3)).symbolSize(0.6*g.x(3)),
	widget.button(stepButton).update(oneStep).size(g.x(3)).symbolSize(0.6*g.x(3))
];

var patternRadio = {id: "patterns", name: "patterns", choices: pattern_list.map(function(d){return d.name}), value:0};
var displayRadio = {id: "display", name: "display", choices: ["world", "weighted sum", "growth", "kernel"], value:0};
//"world", "potential", "growth", "kernel"
var radioWidgets = [
    widget.radio(patternRadio).size(radioBlock.h()).update(selectPattern)
];
var radioWidgets2 = [
    widget.radio(displayRadio).size(radioBlock2.h()).update(switchDisplay).shape("round")
];

var colorToggle = {id:"color", name: "Orli's switch",  value: false};
var hiresToggle = {id:"hires", name: "hi-res",  value: false};
var toggleWidgets = [
	widget.toggle(hiresToggle).update(switchRes).label("bottom").size(10),
	widget.toggle(colorToggle).update(switchColor).label("bottom").size(10)
];

var slider_m = {id:"m", name: "mu", range: [0.01, 0.7], value: 0};
var slider_s = {id:"s", name: "sigma", range: [0.001, 0.07], value: 0};
var slider_R = {id:"R", name: "space", range: [1, 30], value: 10};
var slider_T = {id:"T", name: "time", range: [1, 30], value: 10};
var slider_P = {id:"P", name: "", range: [1, 30], value: 30};
ww = sliderBlock.w();
// width, handleSize, showValue, trackSize, trackBorder, label, fontSize, parameter, name, id, range, value, update, click, X
var sliderWidget_m = widget.slider(slider_m).width(ww).trackSize(8).handleSize(10).update(setParams).showValue(true);
var sliderWidget_s = widget.slider(slider_s).width(ww).trackSize(8).handleSize(10).update(setParams).showValue(true);
var sliderWidget_R = widget.slider(slider_R).width(ww).trackSize(8).handleSize(10).update(setZoom); //.showValue(true);
var sliderWidget_T = widget.slider(slider_T).width(ww).trackSize(8).handleSize(10).update(setParams); //.showValue(true);
var sliderWidget_P = widget.slider(slider_P).width(ww).trackSize(8).handleSize(10).update(setParams);
var sliderWidgets1 = [sliderWidget_s, sliderWidget_m];
var sliderWidgets2 = [sliderWidget_P, sliderWidget_T, sliderWidget_R];

// control setup

var pb = controls.selectAll(".button .play").data(playWidget).enter().append(widget.buttonElement)
	.attr("transform",function(d,i){return "translate("+playBlock.x(0)+","+playBlock.y(0)+")"});	

var bu = controls.selectAll(".button .others").data(buttonWidgets).enter().append(widget.buttonElement)
	.attr("transform",function(d,i){return "translate("+buttonBlock.x(i)+","+buttonBlock.y(0)+")"});	

var rad = controls.selectAll(".radio").data(radioWidgets).enter().append(widget.radioElement)
	.attr("transform",function(d,i){return "translate("+radioBlock.x(0)+","+radioBlock.y(0)+")"});	

var rad2 = controls.selectAll(".radio .two").data(radioWidgets2).enter().append(widget.radioElement)
	.attr("transform",function(d,i){return "translate("+radioBlock2.x(0)+","+radioBlock2.y(0)+")"});	

var tg = controls.selectAll(".toggle").data(toggleWidgets).enter().append(widget.toggleElement)
	.attr("transform",function(d,i){return "translate("+toggleBlock.x(0)+","+toggleBlock.y(i)+")"});	

var sl = controls.selectAll(".slider").data(sliderWidgets1).enter().append(widget.sliderElement)
	.attr("transform",function(d,i){return "translate("+sliderBlock.x(0)+","+sliderBlock.y(i)+")"});

var sl2 = controls.selectAll(".slider .two").data(sliderWidgets2).enter().append(widget.sliderElement)
	.attr("transform",function(d,i){return "translate("+sliderBlock.x(0)+","+sliderBlock.y(i+2)+")"});

sl2.append("text").attr("id","statesLabel").style("font-size",12).attr("transform","translate(0,-16.5)");

//var legend = controls.selectAll(".legend").data([0,0,0]).enter().append("text")
//	.attr("transform",function(d,i){return "translate("+legendBlock.x(0)+","+legendBlock.y(i)+")"})
//  .style("font-size",13);

var colors1 = controls.append("g").attr("class","bars")
	.attr("transform",function(d,i){return "translate("+plotBlock.x(0)+","+plotBlock.y(0)+")"});

var colors2 = controls.append("g").attr("class","bars")
	.attr("transform",function(d,i){return "translate("+plotBlock2.x(0)+","+plotBlock2.y(0)+")"});

var colors3 = controls.append("g").attr("class","bars")
	.attr("transform",function(d,i){return "translate("+plotBlock3.x(0)+","+plotBlock3.y(0)+")"});

var bar1 = colors1.selectAll(".bars").data(d3.range(plotBlock.w()), function(d){return d})
    .enter().append("rect").attr("class", "bars").style("fill", function(d,i){return colormap(d/plotBlock.w())})
    .attr("x", function(d,i){return i}).attr("y", 0).attr("height", 10).attr("width", 1);

var bar2 = colors2.selectAll(".bars").data(d3.range(plotBlock.w()), function(d){return d})
    .enter().append("rect").attr("class", "bars")
    .attr("x", function(d,i){return i}).attr("y", 0).attr("height", 0).attr("width", 1);

var bar3 = colors3.selectAll(".bars").data(d3.range(plotBlock.w()), function(d){return d})
    .enter().append("rect").attr("class", "bars")
    .attr("x", function(d,i){return i}).attr("y", 0).attr("height", 0).attr("width", 1);

var lh = colors1.selectAll(".plot").data(["low","high"]).enter().append("text")
	.text(function(d){return d}).attr("transform", function(d,i){return "translate("+(i*plotBlock.w())+",25)"})
    .style("font-size",12).style("text-anchor","middle");

// maths

const PRECISION = 1000000;
const EPSILON = 1 / PRECISION;
function round(x) { return Math.round(x * PRECISION) / PRECISION; }
function mod(x, n) { return ((x % n) + n) % n; }

// FFT algorithm by Paul Bourke
// http://paulbourke.net/miscellaneous/dft/
function fft1d(dir, re1, im1) {
	/* Do the bit reversal */
	var S = re1.length, m = round(Math.log2(S)), S2 = S >> 1, j1 = 0;
	for (var j=0; j<S-1; j++) {
		if (j < j1) {
			var tmp = re1[j]; re1[j] = re1[j1]; re1[j1] = tmp;
			tmp = im1[j]; im1[j] = im1[j1]; im1[j1] = tmp;
		}
		var j2 = S2;
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
			for (var j=i; j<S; j+=l2) {
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
		var scale_f = 1.0 / S;		
		for (var j=0; j<S; j++) {
			re1[j] *= scale_f;
			im1[j] *= scale_f;
		}
	}
}

function fft2d(dir, re2, im2) {
	var n = re2.length;
    // FFT in 1st direction
	for (var y=0; y<n; y++)
		fft1d(dir, re2[y], im2[y]);
    // transpose
    for (var y=0; y<n; y++)
        for (var x=0; x<y; x++) {
    		[ re2[y][x], re2[x][y] ] = [ re2[x][y], re2[y][x] ];
    		[ im2[y][x], im2[x][y] ] = [ im2[x][y], im2[y][x] ];
        }
    // FFT in 2nd direction
	for (var y=0; y<n; y++)
		fft1d(dir, re2[y], im2[y]);
}

function matrixMult(ar, ai, br, bi, cr, ci) {
	var n = ar.length;
	for (var y=0; y<n; y++) {
		var ar_i = ar[y], ai_i = ai[y];
		var br_i = br[y], bi_i = bi[y];
		var cr_i = cr[y], ci_i = ci[y];
		for (var x=0; x<n; x++) {
			var a = ar_i[x]; var b = ai_i[x];
			var c = br_i[x]; var d = bi_i[x];
			var t = a * (c + d);
			cr_i[x] = t - d*(a+b);
			ci_i[x] = t + c*(b-a);
		}
	}
}

// CA functions

const growthCoreExp  = function(n, m, s, r){return Math.exp( -r*r / (2*s*s) ) * 2 - 1};
const growthCoreStep = function(n, m, s, r){return Math.abs(r-s)<=0.001 ? 0 : r<s ? 1 : -1};

const kernelCoreExp  = function(r){return Math.exp(4 - 1/r/(1-r))};
const kernelCoreStep = function(r){return (r>=1/4 && r<=3/4) ? 1 : 0};

var growthCore = growthCoreExp;
var kernelCore = kernelCoreExp;

function growthFunc(n, m, s) {
	return growthCore(n, m, s, Math.abs(n-m));
}

function kernelFunc(r, b) {
	if (r>=1) return 0;
	var R = r * b.length;
	return kernelCore(R % 1) * b[Math.floor(R)];
}

function calcKernel() {
    var scale = global_R * (isDiscrete ? 1 : worldZoom);
    for (var k=0; k<KN; k++) {
    	var sum = 0.0;
    	var kernel = pattern.kernels[k];
    	for (var y=0; y<N; y++) {
    		for (var x=0; x<N; x++) {
    			var yy = ((y + N/2) % N) - N/2;
    			var xx = ((x + N/2) % N) - N/2;
    			var r = Math.sqrt(xx*xx + yy*yy) / scale / kernel.r;
    			var v = kernelFunc(r, kernel.b);
    			sum += v;
    			Kr[k][y][x] = v;
                // retain kernel for display
    			yy = N - ((y + N/2) % N) - 1;
    			xx = ((x + N/2) % N);
    			K[k][yy][xx] = v;
    		}
            Ki[k][y].fill(0);
    	}

        // pre-calculate F(K)
    	fft2d(1, Kr[k], Ki[k]);

        // normalize kernel
        kernel.sum = sum;
        if (sum > EPSILON) {
        	for (var y=0; y<N; y++) {
        		for (var x=0; x<N; x++) {
        			Kr[k][y][x] /= sum;
        			Ki[k][y][x] /= sum;
        		}
        	}
        }
    }
}

function update(isUpdate) {
    // copy world
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

	// calculate potential U = K * A = F-1( F(K) dot F(A) )
    for (var c=0; c<CN; c++)
        fft2d(1, Ar[c], Ai[c]);
    for (var k=0; k<KN; k++) {
        var kernel = pattern.kernels[k];
        matrixMult(Ar[kernel.c0], Ai[kernel.c0], Kr[k], Ki[k], Ur[k], Ui[k]);
    	fft2d(-1, Ur[k], Ui[k]);
    }

    // calculate growth G = G(U)
    for (var k=0; k<KN; k++) {
        var kernel = pattern.kernels[k];
        var m = KN == 1 ? global_m : kernel.m;
        var s = KN == 1 ? global_s : kernel.s;
    	for (var y=0; y<N; y++) {
    		for (var x=0; x<N; x++) {
    			var u = Ur[k][y][x];
    			G[kernel.c1][y][x] += growthFunc(u, m, s) * kernel.h;
    		}
    	}
    }

    // calculate A = A + dt * G
    var dt = 1/global_T * (isDiscrete ? 1 : worldDelay) * pattern.hScale;
    var mass = 0;
    for (var c=0; c<CN; c++) {
    	for (var y=0; y<N; y++) {
    		for (var x=0; x<N; x++) {
                var g = G[c][y][x];
        		var v = Ao[c][y][x] + dt * g;
                // restrict value to number of states (if not unlimited)
                if (global_P != Math.POSITIVE_INFINITY)
                    v = Math.round(v * global_P) / global_P;
        		if (v<0) v = 0; else if (v>1) v = 1;
        		Ar[c][y][x] = isUpdate ? v : Ao[c][y][x];
                mass += v;
            }
        }
    }

    // statistics
    var scale = global_R * (isDiscrete ? 1 : worldZoom);
    mass = mass / scale / scale;
    if (isUpdate && mass > EPSILON) {
    	gen++;
    	time = round(time + round(dt));
    }
    /*
    var format = d3.format(".2f")
    legend.data([dt, time, mass])
        .text(function(d,i){return i==0 ? "step = "+format(d)+" ms" : i==1 ? "time = "+format(d)+" ms" : "mass = "+format(d)+" mg"});
    */
}

// control functions

// Turbo colormap by Anton Mikhailov
// https://github.com/d3/d3-scale-chromatic/blob/master/src/sequential-multi/turbo.js
function interpolateTurbo(t) {
    t = Math.max(0, Math.min(1, t));
    return "rgb("
        + Math.max(0, Math.min(255, Math.round(34.61 + t * (1172.33 - t * (10793.56 - t * (33300.12 - t * (38394.49 - t * 14825.05))))))) + ", "
        + Math.max(0, Math.min(255, Math.round(23.31 + t * (557.33 + t * (1225.33 - t * (3574.96 - t * (1073.77 + t * 707.56))))))) + ", "
        + Math.max(0, Math.min(255, Math.round(27.2 + t * (3211.1 - t * (15327.97 - t * (27814 - t * (22569.18 - t * 6838.66)))))))
        + ")";
}

function switchRes(d) {
	if (!hiresToggle.value) {
		initWorld(grid=64, delay=0.5, zoom=1);
	} else {
		//initWorld(grid=256, delay=1, zoom=4);
		initWorld(grid=128, delay=1, zoom=2);
	}
    switchDisplay();  // need to re-attach displayArray after switched arrays
	initPattern();
}

function switchColor(d) {
	if (!colorToggle.value) {
		colormap = isDiscrete ? d3.interpolateViridis : d3.interpolateRdYlGn;
        colormapRGB = [
            d3.interpolateRgb.gamma(0.8)("black","rgb(100%,0%,0%)"), 
            d3.interpolateRgb.gamma(0.8)("black","rgb(0%,100%,0%)"), 
            d3.interpolateRgb.gamma(0.8)("black","rgb(0%,0%,100%)")
        ];
	} else {
		//colormap = d3.interpolateRainbow;
        colormap = interpolateTurbo;
        colormapRGB = [
            d3.interpolateRgb.gamma(0.7)("black","rgb(75%,38%,0%)"), 
            d3.interpolateRgb.gamma(0.7)("black","rgb(0%,75%,38%)"), 
            d3.interpolateRgb.gamma(0.7)("black","rgb(38%,0%,75%)")
        ]; // OTV = Orange/Turquoise/Violet
	}

    // draw color bars
    if (CN == 1) {
    	colors1.selectAll(".bars").attr("height", 10).style("fill", function(d,i){return colormap(d/plotBlock.w())});
        colors2.selectAll(".bars").attr("height", 0);
        colors3.selectAll(".bars").attr("height", 0);
    } else {
    	colors1.selectAll(".bars").attr("height", 10).style("fill", function(d,i){return colormapRGB[2](d/plotBlock.w())});
    	colors2.selectAll(".bars").attr("height", 10).style("fill", function(d,i){return colormapRGB[1](d/plotBlock.w())});
    	colors3.selectAll(".bars").attr("height", 10).style("fill", function(d,i){return colormapRGB[0](d/plotBlock.w())});
    }
	draw();
}

function switchDisplay(d) {
	switch (displayRadio.value) {
		case 0: displayArray = Ar; displayN = "CN"; displayMin = 0;  displayMax = 1; break;
		case 1: displayArray = Ur; displayN = "KN"; displayMin = 0;  displayMax = global_m*2; break;
		case 2: displayArray = G;  displayN = "CN"; displayMin = -1; displayMax = 1; break;
		case 3: displayArray = K;  displayN = "KN"; displayMin = 0;  displayMax = 1; break;
	}
	draw();
}

function setParams(d) {
    // get param values from sliders for use in program
    if (KN == 1) {
    	global_m = slider_m.value;
    	if (displayRadio.value == 1) displayMax = global_m * 2;
    	global_s = slider_s.value;
    }
	global_T = slider_T.value;
	global_P = Math.floor(slider_P.value);
    // rightmost of states slider means unlimited 
    if (global_P == slider_P.range[1]) {
        global_P = Math.POSITIVE_INFINITY;
		controls.select("#statesLabel").text(function(d){return "states = unlimited"});
	} else {
		controls.select("#statesLabel").text(function(d){return "states = "+global_P});
	}
    if (playButton.value == 0) {
        update(false);
        draw();
    }
}

function setZoom(d) {
    // re-calc kernel if space changed
	global_R = slider_R.value;
    calcKernel();
    if (playButton.value == 0) {
        update(false);
    	draw();
    }
}

function run() {
    update(true);
    draw();
}

function runPause(d) {
    playButton.value == 1 ? timer = d3.timer(run, 1) : timer.stop();
}

function oneStep(d) {
	if (playButton.value == 1)
		playWidget[0].click();
	run();
}

function selectPattern(d) {
	initPattern(patternRadio.value);
}

function showSliders(isShow){
	if (isShow){
		sl.transition().style("opacity", function(d,i){return 1})
            .select(".track-overlay")
            .style("pointer-events","all");
	} else {
		sl.transition().style("opacity", function(d,i){return i>3 ? 1 : 0})
            .select(".track-overlay")
            .style("pointer-events", function(d,i){return i>3 ? "all" : "none"});
	}
}

// canvas functions

function draw() {
	//canvas.fillStyle = 0;
	//canvas.fillRect(0, 0, worldSize, worldSize);

    // display one channel only
    if (KN == 1) {
    	for (var y=0; y<N; y++) {
            for (var x=0; x<N; x++) {
        		var v = displayArray[0][y][x];
        		v = (v - displayMin) / (displayMax - displayMin);
                canvas.fillStyle = colormap(v);
                canvas.fillRect(xScale(x), yScale(y), pixelSize, pixelSize);
            }
        }
    // display for multi channels
    } else if (displayN == "CN") {
    	for (var y=0; y<N; y++) {
            for (var x=0; x<N; x++) {
                var color = d3.rgb(0, 255*0.2, 0);
                // accumulate color for each channel
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
    // display for multi kernels
    } else if (displayN == "KN") {
    	for (var y=0; y<N; y++) {
            for (var x=0; x<N; x++) {
                var color = d3.rgb(0, 255*0.2, 0);
                // accumulate color for each kernel
            	for (var k=0; k<KN; k++) {
                    var kernel = pattern.kernels[k];
            		var v = displayArray[k][y][x] / 2;  // scale value for better display
                    color0 = d3.rgb(colormapRGB[kernel.c0](v));  // use source channel colormap
                    color.r += color0.r;
                    color.g += color0.g;
                    color.b += color0.b;
                }
                canvas.fillStyle = color;
                canvas.fillRect(xScale(x), yScale(y), pixelSize, pixelSize);
            }
        }
    }

    if (displayRadio.value == 3)
        drawLines();
}

function drawLines() {
    // draw lines of growth function
    var panelH = 30;
    var panelX = 0;
    var panelY = panelH*2+5;
    canvas.fillStyle = "rgba(255,255,255,0.5)";
    canvas.fillRect(0, 0, worldSize, panelH+5);
    canvas.fillRect(0, panelH+6, worldSize, panelH+5);
	for (var k=0; k<KN; k++) {
        canvas.beginPath();
        canvas.moveTo(-10, panelY);
        var kernel = pattern.kernels[k];
        for (var x=0; x<=worldSize; x++) {
            var m = KN == 1 ? global_m : kernel.m;
            var s = KN == 1 ? global_s : kernel.s;
            var v = growthFunc(x/worldSize, m, s) * kernel.h;
            canvas.lineTo(panelX + x, panelY - (v+1)*panelH);
        }
        canvas.lineWidth = 1;
        canvas.strokeStyle = KN == 1 ? "black" : colormapRGB[kernel.c0](1);
        canvas.stroke();
    }

    // draw lines of kernel function
    panelH = 40;
    panelX = worldSize/2 + pixelSize/2;
    panelY = worldSize-30 - pixelSize/2;
    var scale = global_R * (isDiscrete ? 1 : worldZoom);
	for (var k=0; k<KN; k++) {
        canvas.beginPath();
        canvas.moveTo(-10, panelY - pixelSize/2);
        var kernel = pattern.kernels[k];
        for (var x=-worldSize; x<=worldSize; x++) {
            var v = kernelFunc(Math.abs(x / worldSize * N/2 / scale / kernel.r), kernel.b);
            canvas.lineTo(panelX + x/2, panelY - v*panelH);
        }
        canvas.lineWidth = 2;
        var color = d3.rgb(KN == 1 ? "white" : colormapRGB[kernel.c0](1));
        color.opacity = 0.7;
        canvas.strokeStyle = color;
        canvas.stroke();
    }
}

// world & pattern functions

function clearWorld() {
	for (var c=0; c<CN; c++)
        for (var y=0; y<N; y++) {
            Ar[c][y].fill(0);
            Ai[c][y].fill(0);
        }
}

function placePattern(p, shifty, shiftx, angle, zoom) {
    var ch = p.cells[0].length;
    var cw = p.cells[0][0].length;
    // scale pattern
    var scale = p.displayR / p.R * zoom * (isDiscrete ? 1 : worldZoom);
    // rotate pattern
	var sin = Math.sin(angle / 180 * Math.PI);
	var cos = Math.cos(angle / 180 * Math.PI);
	var h = (Math.abs(ch*cos) + Math.abs(cw*sin) + 1) * scale - 1;
	var w = (Math.abs(cw*cos) + Math.abs(ch*sin) + 1) * scale - 1;
    // shift pattern
    var y0 = Math.floor(N/2 - h/2 + shifty);
    var x0 = Math.floor(N/2 - w/2 + shiftx);
    // copy intrapolated cells from pattern to world
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

function initPattern(id, isInit=false) {
	if (id>=0) {
		pattern = pattern_list[id];
        CN = pattern.cells.length;
        KN = pattern.kernels.length;
    }

    initParams();
    initCells(isInit);
}

function initParams() {
    // use step kernel for discrete CA e.g. Conway's Game of Life
    isDiscrete = pattern.R<=2;
    growthCore = growthCoreExp;  // isDiscrete ? growthCoreStep : growthCoreExp;
    kernelCore = isDiscrete ? kernelCoreStep : kernelCoreExp;
    switchColor();

    if (KN == 1) {
	    var kernel = pattern.kernels[0];
    	sliderWidget_m.click(kernel.m);
    	sliderWidget_s.click(kernel.s);
        showSliders(true);
    } else {
        // hide mu & sigma when too many kernels 
        showSliders(false);
    }

	sliderWidget_T.click(pattern.T);
	sliderWidget_P.click(isDiscrete ? 1 : slider_P.range[1]);
}

function initCells(isInit=false) {
    clearWorld();
	if (isInit==true || isDiscrete) {
        // place 1 pattern at center
        var zoom = 1;
        placePattern(pattern, 0, 0, 0, zoom);
        sliderWidget_R.click(zoom*pattern.displayR);
    } else {
        // random place 1 to 3 patterns
        randomPlacePattern(pattern);
    }
    gen = 0;
    time = 0.0;
	update(false);
    draw();
}

function randomPlacePattern(p) {
	var chance = Math.random();
    var zoom = 1;
	if (p.name.includes("(s)") && chance>0.7) {
        // random place 3 patterns
		zoom = 0.7;
		placePattern(p, -N/4, 0, Math.random()*360, zoom);
		placePattern(p, +N/4, +N/4, Math.random()*360, zoom);
		placePattern(p, +N/4, -N/4, Math.random()*360, zoom);
	} else if (p.name.includes("(s)") && chance>0.3) {
        // random place 2 patterns
		zoom = 0.8;
		placePattern(p, +N/4, +N/4, Math.random()*360, zoom);
		placePattern(p, -N/4, -N/4, Math.random()*360, zoom);
	} else {
        // random place 1 pattern
		placePattern(p, 0, 0, Math.random()*360, zoom);
	}

    sliderWidget_R.click(zoom*pattern.displayR);
}

// init functions

function initWorld(grid, delay, zoom) {
	N = grid;
	pixelSize = Math.ceil(worldSize / N);
	worldDelay = delay;
	worldZoom = zoom;

    // choose array for grid size
    Ar = setAr[N];
    Ai = setAi[N];
    Ao = setAo[N];
    K  = setK [N];
    Kr = setKr[N];
    Ki = setKi[N];
    Ur = setUr[N];
    Ui = setUi[N];
    G  = setG [N];
	displayArray = Ar;
	displayMin = 0;
	displayMax = 1;
    displayN = "CN";

    // grid to canvas conversion tool
	xScale = d3.scaleLinear().domain([0,N]).range([0,worldSize]);
	yScale = d3.scaleLinear().domain([0,N]).range([0,worldSize]);
}

function initAllArrays() {
	setAr = initArraySet(3);
	setAi = initArraySet(3);
	setAo = initArraySet(3);
	setK  = initArraySet(16);
	setKr = initArraySet(16);
	setKi = initArraySet(16);
	setUr = initArraySet(16);
	setUi = initArraySet(16);
	setG  = initArraySet(3);
}

function initArraySet(channels) {
    a = [];
    // init array for each grid size
    a[64] = initArray(64, channels);
    a[128] = initArray(128, channels);
    return a;
}

function initArray(size, channels) {
    var arr = [];
	for (var j=0; j<channels; j++) {
        var ch = [];
        for (var i=0; i<size; i++)
            ch.push(new Array(size).fill(0));
        arr.push(ch);
    }
	return arr;
}

initAllArrays();
initWorld(grid=64, delay=0.5, zoom=1);
initPattern(id=0, isInit=true);
//console.log(Object.getOwnPropertyNames(sliderWidget_R));

})()
