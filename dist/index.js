"use strict";
// var addon = require("./build/Release/napi-addon-fdtd.node");
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
var bindings_1 = __importDefault(require("bindings"));
var addon = (0, bindings_1.default)('napi-addon-fdtd');
exports.default = addon;
