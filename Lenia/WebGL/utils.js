//from http://webglfundamentals.org/webgl/lessons/webgl-boilerplate.html
"use strict";

function compileShader(gl, shaderSource, shaderType) {
    var shader = gl.createShader(shaderType);
    gl.shaderSource(shader, shaderSource);
    gl.compileShader(shader);
    var success = gl.getShaderParameter(shader, gl.COMPILE_STATUS);
    if (!success) {
        throw "failed to compile shader: " + gl.getShaderInfoLog(shader);
    }
    return shader;
}

function createProgramFromSources(gl, vertexShaderSource, fragmentShaderSource) {
    var vertexShader = compileShader(gl, vertexShaderSource, gl.VERTEX_SHADER);
    var fragmentShader = compileShader(gl, fragmentShaderSource, gl.FRAGMENT_SHADER);

    var program = gl.createProgram();

    gl.attachShader(program, vertexShader);
    gl.attachShader(program, fragmentShader);
    gl.linkProgram(program);

    var success = gl.getProgramParameter(program, gl.LINK_STATUS);
    if (!success) {
        throw "filed to link program: " + gl.getProgramInfoLog (program);
    }

    return program;
}



// from https://stackoverflow.com/questions/4878145/javascript-and-webgl-external-scripts

function loadFile(url, urlIndex, callback) {
    var request = new XMLHttpRequest();
    request.open('GET', url, true);
    //request.setRequestHeader("Accept", "application/vnd.github.3.raw");

    request.onreadystatechange = function () {
        if (request.readyState == 4) {
            if (request.status == 0 || request.status == 200) {
                callback(request.responseText, urlIndex)
            } else {
                throw "failed to download shader: (" + url + ") " + request.statusText;
            }
        }
    };

    request.send(null);    
}

function loadShaderFiles(urls, callback) {
    var numUrls = urls.length;
    var numComplete = 0;
    var results = [];

    function partialCallback(text, urlIndex) {
        results[urlIndex] = text;
        numComplete++;

        if (numComplete == numUrls) {
            callback(results);
        }
    }

    for (var i = 0; i < numUrls; i++) {
        loadFile(urls[i], i, partialCallback);
    }
}
