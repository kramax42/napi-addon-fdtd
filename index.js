var addon = require("./build/Release/napi-addon-fdtd.node");


module.exports = addon;


let data = addon.getData2D([1, 10, 1], true, [1, 1.2], 2);

for (let j = 0; j < 50; ++j) {
  data = addon.getData2D([1, 10, 1], false, [1, 1.2], 2);

}

const fs = require("fs");
const path = require("path");
fs.writeFileSync(path.resolve(__dirname, "tmp.txt"), JSON.stringify(data.dataHy), function (err) {
  if (err) {
    return console.log(err);
  }
  console.log("The file was saved!");
}); // Orfs.writeFileSync('/tmp/test-sync', 'Hey there!');


// let data = addon.getData3D([1, 10], true, [1,2,1,1], 2, 0);
// for (let j = 0; j < 150; ++j) {
//   data = addon.getData3D([1, 10], false, [1,2,1,1], 2, 0);
// }

// const fs = require("fs");
// const path = require("path");
// fs.writeFileSync(path.resolve(__dirname, "tmp.txt"), JSON.stringify(data.dataY), function (err) {
//   if (err) {
//     return console.log(err);
//   }
//   console.log("The file was saved!");
// }); // Orfs.writeFileSync('/tmp/test-sync', 'Hey there!');