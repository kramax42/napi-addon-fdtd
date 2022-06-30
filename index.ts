// Import addon without 'bindings' package
// let addon = require("./build/Release/napi-addon-fdtd.node");

import bindings from 'bindings';
const addon = bindings('napi-addon-fdtd')

export default addon;


