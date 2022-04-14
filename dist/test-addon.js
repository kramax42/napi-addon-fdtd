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
    var sigma = [0, 0.04];
    var srcPosition = [0.4, 0.8];
    var data = index_1.default.getData2D(condition, true, eps, 2, srcPosition, sigma);
    for (var j = 0; j < 50; ++j) {
        data = index_1.default.getData2D(condition, false, eps, 2, srcPosition, sigma);
    }
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
    var eps = Array(epsSize * 2).fill(0).map(function (_) { return Math.random() * 10; });
    var condition = [1, 10];
    var reload = true;
    var data = index_1.default.getData3D(condition, reload, [4, 5, 6, 7], 2, 0);
    reload = false;
    for (var j = 0; j < 150; ++j) {
        //eps, epsSize
        data = index_1.default.getData3D(condition, reload, [4, 5, 6, 7], 2, 0);
    }
    fs_1.default.writeFileSync(path_1.default.resolve(__dirname, "tmp.txt"), JSON.stringify(data.dataY), 
    // @ts-ignore
    function (err) {
        if (err) {
            return console.log(err);
        }
        console.log("The file was saved!");
    }); // Orfs.writeFileSync('/tmp/test-sync', 'Hey there!');
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
