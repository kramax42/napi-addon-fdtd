
type GetData2D = (
  condition: [number, number],
  reload: boolean,
  materialMatrix: number[],
  rows: number,
  eps: number[],
  mu: number[],
  sigma: number[],
  dataToReturn: number,
  srcPositionRelativeSet: number[]
) => {
  dataX: number[],
  dataY: number[],
  dataEz: number[],
  rows: number,
  cols: number,
  timeStep: number,
  max: number,
  min: number,
}

type GetData1D = (
  condition: number[],
  reload: boolean,
  eps: number[],
  epsSize: number,
  sourcePositionRelative: number[],
  sigma: number[]
  // dataReturnType: number,
) => any

// interface Module {
//     sayHi: (id: number) => void;
//     getData3D: GetData3D;
//     getData2D: GetData2D;
//   }
  
  // declare module 'bindings' {
  //   export default function(string: 'hello'): Module;
  // }

  declare module 'napi-addon-fdtd' {
    export const getData2D: GetData2D;
    export const getData1D: GetData1D;
  }