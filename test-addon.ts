import addon from "./index";

const test1D = () => {
  const omega = 3.5;
  const tau = 8;
  const materialVector = [1, 0, 2, 0, 1];
  const eps = [1.0, 1.2, 1.1];
  const mu = [0.51, 0.5, 0.57];
  const sigma = [1.0, 0.001, 1.0];
  const srcPosition = [0.5];

  const isReload = true;

  let fdtd = new addon.Fdtd1D({
    omega,
    tau,
    isReload,
    materialVector,
    eps,
    mu,
    sigma,
    srcPosition
  });


  let data;
  for (let j = 0; j < 49; ++j) {
    data = fdtd.getNextTimeLayer();
  }
  console.log(data);
};

const test2D = () => {
  const lambda = 1;
  const beamsize = 1;
  const materialVector = [1, 0, 2, 0];
  const eps = [1.0, 1.2, 1.1];
  const mu = [0.51, 0.5, 0.57];
  const sigma = [1.0, 0.001, 1.0];
  const srcPosition = [0.1, 0.1];
  const dataReturnType = 0;
  let isReload = true;

  let fdtd = new addon.Fdtd2D({
    lambda,
    beamsize,
    isReload,
    materialVector,
    eps,
    mu,
    sigma,
    dataReturnType,
    srcPosition
  });

  let data;
  for (let j = 0; j < 49; ++j) {
    data = fdtd.getNextTimeLayer();
  }
  console.log(data);
};

function testMemoryUsage() {
  const used = process.memoryUsage().heapUsed / 1024 / 1024;
  console.log(
    `The script uses approximately ${Math.round(used * 100) / 100} MB`
  );
}


test1D();
// test2D();
// testMemoryUsage();
