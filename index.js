var addon = require("./build/Release/napi-addon-fdtd.node");

module.exports = addon;

const fs = require("fs");
const path = require("path");

const test1D = () => {
  const condition = [1, 10, 1];
  const eps = [1, 1.2];
  const sigma = [0, 0.04];
  const srcPosition = [0.4, 0.8];

  let data = addon.getData2D(condition, true, eps, 2, srcPosition, sigma);

  for (let j = 0; j < 150; ++j) {
    data = addon.getData2D(condition, false, eps, 2, srcPosition, sigma);
  }

  fs.writeFileSync(
    path.resolve(__dirname, "tmp.txt"),
    JSON.stringify(data.dataHy),
    function (err) {
      if (err) {
        return console.log(err);
      }
      console.log("The file was saved!");
    }
  ); // Orfs.writeFileSync('/tmp/test-sync', 'Hey there!');
};

const test2D = () => {
  let data = addon.getData3D([1, 10], true, [1, 2, 1, 1], 2, 0);
  for (let j = 0; j < 150; ++j) {
    data = addon.getData3D([1, 10], false, [1, 2, 1, 1], 2, 0);
  }

  fs.writeFileSync(
    path.resolve(__dirname, "tmp.txt"),
    JSON.stringify(data.dataY),
    function (err) {
      if (err) {
        return console.log(err);
      }
      console.log("The file was saved!");
    }
  ); // Orfs.writeFileSync('/tmp/test-sync', 'Hey there!');
};


test1D();
test2D();