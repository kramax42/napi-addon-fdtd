import addon from "./index";

const test1D = () => {
  const condition = [3.5, 8];
  const materialVector = [1, 0, 2, 0, 1];
  const eps = [1.0, 1.2, 1.1];
  const mu = [0.51, 0.5, 0.57];
  const sigma = [1.0, 0.001, 1.0];
  const srcPosition = [0.5];

  let reload = true;
  let data = addon.getData1D(
    condition,
    reload,
    materialVector,
    materialVector.length,
    eps,
    mu,
    sigma,
    srcPosition
  );

  reload = false;
  for (let j = 0; j < 5; ++j) {
    data = addon.getData1D(
      condition,
      reload,
      materialVector,
      materialVector.length,
      eps,
      mu,
      sigma,
      srcPosition
    );
  }
  console.log(data);
};

const test2D = () => {
  const condition = [1, 10];
  const materialMatrix = [1, 0, 2, 0];
  const matrixSize = 2;
  const eps = [1.0, 1.2, 1.1];
  const mu = [0.51, 0.5, 0.57];
  const sigma = [1.0, 0.001, 1.0];
  const returnDataNumber = 0;
  const srcPosition = [0.1, 0.1];

  let reload = true;
  let data = addon.getData2D(
    condition,
    reload,
    materialMatrix,
    matrixSize,
    eps,
    mu,
    sigma,
    returnDataNumber,
    srcPosition
  );

  reload = false;
  for (let j = 0; j < 50; ++j) {
    data = addon.getData2D(
      condition,
      reload,
      materialMatrix,
      matrixSize,
      eps,
      mu,
      sigma,
      returnDataNumber,
      srcPosition
    );
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
test2D();
// testMemoryUsage();
