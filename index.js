var addon = require("./build/Release/napi-addon-fdtd.node");

module.exports = addon;

// let data = addon.getFdtd2D([1, 10, 1], true);
// for (let j = 0; j < 150; ++j) {
//   data = addon.getFdtd2D([1, 10, 1], false);
// }

// const fs = require("fs");
// const path = require("path");
// fs.writeFileSync(path.resolve(__dirname, "tmp.txt"), JSON.stringify(data.dataY), function (err) {
//   if (err) {
//     return console.log(err);
//   }
//   console.log("The file was saved!");
// }); // Orfs.writeFileSync('/tmp/test-sync', 'Hey there!');
