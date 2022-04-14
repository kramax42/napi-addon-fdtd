
type GetData3D = (
    condition: number[],
    reload: boolean,
    eps: number[],
    epsSize: number,
    dataReturnType: number,
) => any

type GetData2D = (
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
    export const getData3D: GetData3D;
    export const getData2D: GetData2D;
  }