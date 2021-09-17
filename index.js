//const addon = require("bindings")("FDTD.node");
const addon = require("./build/Release/napi-addon-fdtd.node");

// console.log("exports: ", addon);
// exports.addonFDTD = addon;
module.exports = addon;

  // addon.getFDTD_2D([1, 10, 1], false);
  // addon.getFDTD_2D([1, 10, 1], false);
  // addon.getFDTD_2D([1, 10, 1], false);
  //  console.log("data1: ", addon.getFDTD_2D([1, 10, 1], false).currentTick);
  // console.log("data1: ", addon.getFDTD_2D([1, 10, 1], false));
  // console.log("data1: ", addon.getFDTD_2D([1, 10, 1], true));

//  console.log("data2: ", addon.getFDTD_2D([1, 11, 1]));

//  console.log("test: ");
//  console.log("data1: ", addon.getFDTD_3D([1, 3, 1], false));
//  console.log("data2: ", addon.getFDTD_3D([1, 3, 1], false));
//  console.log("data3: ", addon.getFDTD_3D([1, 3, 1], false));
//  console.log("data4: ", addon.getFDTD_3D([1, 3, 1], true));

