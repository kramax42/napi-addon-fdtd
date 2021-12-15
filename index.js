var addon = require("./build/Release/napi-addon-fdtd.node");

module.exports = addon;
// console.log(addon);
// console.log("data1: ", addon.getFDTD_2D([1, 10, 1], false));
console.log("data1: ", addon.getFDTD_2D([1, 10, 1], true));

const conditions = [1, 3, 1, 1.5];
const reload = true;
const refrMatrix = [1, 2, 1, 1];
const refrMatrixSize = 2;
const dataType = 0; // Ez
console.log(
  "data1: ",
  addon.getFDTD_3D_DIFRACTION(
    conditions,
    reload,
    refrMatrix,
    refrMatrixSize,
    dataType
  )
);

for (let j = 0; j < 20; ++j) {
  addon.getFDTD_3D_DIFRACTION(
    conditions,
    false,
    refrMatrix,
    refrMatrixSize,
    dataType
  );
}

let {
  dataX,
  dataY,
  dataEz: dataVal,
} = addon.getFDTD_3D_DIFRACTION(
  conditions,
  reload,
  refrMatrix,
  refrMatrixSize,
  dataType
);

//   console.log("data2: ", addon.getFDTD_3D_DIFRACTION([1, 3, 1, 1.5], false));
//  console.log("data3: ", addon.getFDTD_3D([1, 3, 1], false));
// //  console.log("data4: ", addon.getFDTD_3D([1, 3, 1], true));

//   console.log("data5: ", addon.getFDTD_3D_INTERFERENCE([1, 3, 1], false));
//   console.log("data6: ", addon.getFDTD_3D_INTERFERENCE([1, 3, 1], false));

const { createCanvas } = require("canvas");

const width = 600;
const height = 600;

const canvas = createCanvas(width, height);
const ctx = canvas.getContext("2d");

let brushSize = width / 285;
let brushBlurSize = brushSize * 1.15;

var r = brushSize + brushBlurSize;
var d = r * 2;

const brushCanvas = createCanvas(d, d);
const brushctx = brushCanvas.getContext("2d");

brushctx.shadowOffsetX = d;
brushctx.shadowBlur = brushBlurSize;
brushctx.shadowColor = "black";

// draw circle in the left to the canvas

brushctx.beginPath();
brushctx.arc(-r, r, brushSize, 0, Math.PI * 2, true);
brushctx.closePath();
brushctx.fill();

///////////////////////
const gridSizeFromBackend = 320;

// for (let i = 0; i < dataX.length; ++i) {
//   data.push([
//     (dataX[i] * width) / gridSizeFromBackend,
//     (dataY[i] * height) / gridSizeFromBackend,
//     (dataVal[i] + Math.abs(minVal)) / (maxVal + Math.abs(minVal)),
//   ]);
// }

var data = [];
for (var i = 0; i < 1000; ++i) {
  data.push([Math.random() * width, Math.random() * height, Math.random()]);
}

////////////////////////////////

const minVal = -1;
const maxVal = 1;

for (var i = 0; i < dataX.length; ++i) {
  var x = (dataX[i] * width) / gridSizeFromBackend;
  var y = (dataY[i] * height) / gridSizeFromBackend;
  var alpha = (dataVal[i] + Math.abs(minVal)) / (maxVal + Math.abs(minVal)); // using value as alpha

  // draw with the circle brush with alpha

  ctx.globalAlpha = alpha;
  ctx.drawImage(brushCanvas, x - r, y - r);
}

////////////////// Gradientvar levels = 256;
const levels = 256;
const gradientCanvas = createCanvas(10, levels);

var ctxGradient = gradientCanvas.getContext("2d");

var gradientColors = {
  0.4: "blue",
  0.5: "cyan",
  0.6: "lime",
  0.8: "yellow",
  1.0: "red",
};

// add color to gradient stops

var gradient = ctxGradient.createLinearGradient(0, 0, 0, levels);
for (var pos in gradientColors) {
  gradient.addColorStop(+pos, gradientColors[pos]);
}

ctxGradient.fillStyle = gradient;
ctxGradient.fillRect(0, 0, 1, levels);
var gradientPixels = ctxGradient.getImageData(0, 0, 1, levels).data;

////end gradient

var imageData = ctx.getImageData(0, 0, width, height);
var pixels = imageData.data;
var len = pixels.length / 4;
while (len--) {
  var id = len * 4 + 3;
  var alpha = pixels[id] / 256; // why not `gradientLevels`?

  var colorOffset = Math.floor(alpha * (levels - 1));
  pixels[id - 3] = gradient[colorOffset * 4]; // red

  pixels[id - 2] = gradient[colorOffset * 4 + 1]; // green

  pixels[id - 1] = gradient[colorOffset * 4 + 2]; // blue
}
ctx.putImageData(imageData, 0, 0);

const fs = require("fs");
const buffer = canvas.toBuffer("image/png");
fs.writeFileSync("./test.png", buffer);
// })
