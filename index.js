const addon = require("./build/Release/napi-addon-fdtd.node");

module.exports = addon;
console.log(addon);
// console.log("data1: ", addon.getFDTD_2D([1, 10, 1], false));
// console.log("data1: ", addon.getFDTD_2D([1, 10, 1], true));

  console.log("data1: ", addon.getFDTD_3D_DIFRACTION([1, 3, 1, 1.5], false));
  console.log("data2: ", addon.getFDTD_3D_DIFRACTION([1, 3, 1, 1.5], false));
//  console.log("data3: ", addon.getFDTD_3D([1, 3, 1], false));
//  console.log("data4: ", addon.getFDTD_3D([1, 3, 1], true));

