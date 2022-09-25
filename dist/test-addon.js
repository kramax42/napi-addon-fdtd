"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
var index_1 = __importDefault(require("./index"));
var test1D = function () {
    var condition = [3.5, 8];
    var materialVector = [1, 0, 2, 0, 1];
    var eps = [1.0, 1.2, 1.1];
    var mu = [0.51, 0.5, 0.57];
    var sigma = [1.0, 0.001, 1.0];
    var srcPosition = [0.5];
    var reload = true;
    var fdtd = new index_1.default.Fdtd1D(condition, reload, materialVector, materialVector.length, eps, mu, sigma, srcPosition);
    var data;
    for (var j = 0; j < 49; ++j) {
        data = fdtd.getNextTimeLayer();
    }
    console.log(data);
};
var test2D = function () {
    var condition = [1, 10];
    var materialMatrix = [1, 0, 2, 0];
    var matrixSize = 2;
    var eps = [1.0, 1.2, 1.1];
    var mu = [0.51, 0.5, 0.57];
    var sigma = [1.0, 0.001, 1.0];
    var returnDataNumber = 0;
    var srcPosition = [0.1, 0.1];
    var reload = true;
    var data = index_1.default.getData2D(condition, reload, materialMatrix, matrixSize, eps, mu, sigma, returnDataNumber, srcPosition);
    reload = false;
    for (var j = 0; j < 50; ++j) {
        data = index_1.default.getData2D(condition, reload, materialMatrix, matrixSize, eps, mu, sigma, returnDataNumber, srcPosition);
    }
    console.log(data);
};
function testMemoryUsage() {
    var used = process.memoryUsage().heapUsed / 1024 / 1024;
    console.log("The script uses approximately ".concat(Math.round(used * 100) / 100, " MB"));
}
test1D();
// test2D();
// testMemoryUsage();
