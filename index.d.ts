
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
  condition: [number, number],
  reload: boolean,
  materialVector: number[],
  rows: number,
  eps: number[],
  mu: number[],
  sigma: number[],
  srcPositionRelativeSet: number[]
) => any

// interface Module {
//     sayHi: (id: number) => void;
//     getData3D: GetData3D;
//     getData2D: GetData2D;
//   }
  
  // declare module 'bindings' {
  //   export default function(string: 'hello'): Module;
  // }


  type Fdtd1dOptions = {
    omega: number;
    tau: number;
    isReload: boolean;
    materialVector: number[];
    eps: number[];
    mu: number[];
    sigma: number[];
    srcPosition: number[];
  };

  type Fdtd1dOutput = {
    max: number;
    min: number;
    dataX: number[];
    dataY: number[];
  };

  type Fdtd2dOptions = {
    lambda: number;
    beamsize: number;
    isReload: boolean;
    materialVector: number[];
    eps: number[];
    mu: number[];
    sigma: number[];
    dataReturnType: number;
    srcPosition: number[];
  };

  type Fdtd2dOutput = {
    max: number;
    min: number;
    dataX: number[];
    dataY: number[];
    rows: number;
    cols: number;  
    timestep: number;  
    dataEz?: number[];
    dataHx?: number[];
    dataHy?: number[];
    dataEnergy?: number[];
  };

  declare module 'napi-addon-fdtd' {
    // export const getData2D: GetData2D;
    // export const getData1D: GetData1D;
    class Fdtd1D {
      constructor(options: Fdtd1dOptions);
      getNextTimeLayer(): Fdtd1dOutput;
    }

    class Fdtd2D {
      constructor(options: Fdtd2dOptions);
      getNextTimeLayer(): Fdtd2dOutput;
    }

    class Fdtd2DTFSF {
      constructor();
      getNextTimeLayer(): Fdtd2dOutput;
    }
    
  }


  