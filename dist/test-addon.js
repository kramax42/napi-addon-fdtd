"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
var fs_1 = __importDefault(require("fs"));
var path_1 = __importDefault(require("path"));
var index_1 = __importDefault(require("./index"));
var test1D = function () {
    var condition = [1, 10, 1];
    var eps = [1, 1.2];
    var materialSize = 2;
    var sigma = [0, 0.04];
    var srcPosition = [0.4, 0.8];
    var data = index_1.default.getData1D(condition, true, eps, materialSize, srcPosition, sigma);
    for (var j = 0; j < 150; ++j) {
        data = index_1.default.getData1D(condition, false, eps, materialSize, srcPosition, sigma);
    }
    console.log(data);
    fs_1.default.writeFileSync(path_1.default.resolve(__dirname, "tmp.txt"), JSON.stringify(data.dataHy), 
    // @ts-ignore
    function (err) {
        if (err) {
            return console.log(err);
        }
        console.log("The file was saved!");
    }); // Orfs.writeFileSync('/tmp/test-sync', 'Hey there!');
};
var test2D = function () {
    var epsSize = 40;
    var condition = [1, 10];
    var materialMatrix = [1, 0, 2, 0];
    var matrixSize = 2;
    var eps = [1.0, 1.2, 1.1];
    var mu = [0.51, 0.5, 0.57];
    var sigma = [1.0, 0.001, 1.0];
    var returnDataNumber = 0;
    var srcPosition = [0, 0];
    var reload = true;
    var data = index_1.default.getData2D(condition, reload, materialMatrix, matrixSize, eps, mu, sigma, returnDataNumber, srcPosition);
    reload = false;
    for (var j = 0; j < 50; ++j) {
        data = index_1.default.getData2D(condition, reload, materialMatrix, matrixSize, eps, mu, sigma, returnDataNumber, srcPosition);
    }
    console.log(data);
    // fs.writeFileSync(
    //   path.resolve(__dirname, "tmp.txt"),
    //   JSON.stringify(data.dataY),
    //    // @ts-ignore
    //   function (err) {
    //     if (err) {
    //       return console.log(err);
    //     }
    //     console.log("The file was saved!");
    //   }
    // ); // Orfs.writeFileSync('/tmp/test-sync', 'Hey there!');
};
function testMemoryUsage() {
    // const arr = [1, 2, 3, 4, 5, 6, 9, 7, 8, 9, 10];
    // const arr = Array(1e7).fill(1e3);
    // arr.reverse();
    var used = process.memoryUsage().heapUsed / 1024 / 1024;
    console.log("The script uses approximately ".concat(Math.round(used * 100) / 100, " MB"));
}
test1D();
// test2D();
// testMemoryUsage();
