"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
var index_1 = __importDefault(require("./index"));
var test1D = function () {
    var omega = 3.5;
    var tau = 8;
    var materialVector = [1, 0, 2, 0, 1];
    var eps = [1.0, 1.2, 1.1];
    var mu = [0.51, 0.5, 0.57];
    var sigma = [1.0, 0.001, 1.0];
    var srcPosition = [0.5];
    var isReload = true;
    var fdtd = new index_1.default.Fdtd1D({
        omega: omega,
        tau: tau,
        isReload: isReload,
        materialVector: materialVector,
        eps: eps,
        mu: mu,
        sigma: sigma,
        srcPosition: srcPosition
    });
    var data;
    for (var j = 0; j < 49; ++j) {
        data = fdtd.getNextTimeLayer();
    }
    console.log(data);
};
var test2D = function () {
    var lambda = 1;
    var beamsize = 1;
    // const materialVector = [1, 0, 2, 0, 1, 1, 1,1 ,1];
    // const materialVector = [1, 0, 1, 0, 1, 1, 1, 1 ,1];
    var materialVector = [0, 0, 0, 0];
    // const eps = [1.0, 1.2, 1.1];
    // const mu = [0.51, 0.5, 0.57];
    // const sigma = [1.0, 0.001, 1.0];
    var eps = [1.2, 1.15];
    var mu = [1, 1];
    var sigma = [0, 1e+7];
    var srcPosition = [0.1, 0.1];
    var dataReturnType = 0;
    var isReload = true;
    var fdtd = new index_1.default.Fdtd2D({
        lambda: lambda,
        beamsize: beamsize,
        isReload: isReload,
        materialVector: materialVector,
        eps: eps,
        mu: mu,
        sigma: sigma,
        dataReturnType: dataReturnType,
        srcPosition: srcPosition
    });
    var data;
    for (var j = 0; j < 300; ++j) {
        data = fdtd.getNextTimeLayer();
    }
    console.log(data);
};
function testMemoryUsage() {
    var used = process.memoryUsage().heapUsed / 1024 / 1024;
    console.log("The script uses approximately ".concat(Math.round(used * 100) / 100, " MB"));
}
var test2DTFSF = function () {
    var fdtd = new index_1.default.Fdtd2DTFSF();
    var data;
    for (var j = 0; j < 30; ++j) {
        data = fdtd.getNextTimeLayer();
    }
    // console.log(fdtd);
    console.log(data);
};
// test1D();
// test2D();
test2DTFSF();
// testMemoryUsage();
