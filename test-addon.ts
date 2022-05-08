import fs from 'fs';
import path from 'path';

import addon from './index';



const test1D = () => {

  const condition = [1, 10, 1];
  const eps = [1, 1.2];
  const materialSize = 2;
  const sigma = [0, 0.04];
  const srcPosition = [0.4, 0.8];

  let data = addon.getData1D(condition, true, eps, materialSize, srcPosition, sigma);

  for (let j = 0; j < 150; ++j) {
    data = addon.getData1D(condition, false, eps, materialSize, srcPosition, sigma);
  }

  console.log(data)

  fs.writeFileSync(
    path.resolve(__dirname, "tmp.txt"),
    JSON.stringify(data.dataHy),
    // @ts-ignore
    function (err) {
      if (err) {
        return console.log(err);
      }
      console.log("The file was saved!");
    }
  ); // Orfs.writeFileSync('/tmp/test-sync', 'Hey there!');
};

const test2D = () => {
  const epsSize = 40;

  const condition = [1, 10]
  const materialMatrix = [1,0,2,0]
  const matrixSize = 2;
  const eps = [1.0, 1.2, 1.1];
  const mu = [0.51, 0.5, 0.57];
  const sigma = [1.0, 0.001, 1.0];
  const returnDataNumber = 0;
  const srcPosition = [0,0];

  let reload = true;
  let data = addon.getData2D(condition, reload, materialMatrix, matrixSize, eps, mu, sigma, returnDataNumber, srcPosition);

  reload = false;
  for (let j = 0; j < 200; ++j) {
    data = addon.getData2D(condition, reload, materialMatrix, matrixSize, eps, mu, sigma, returnDataNumber, srcPosition);
  }

  console.log(data)
  
  // fs.writeFileSync(
  //   path.resolve(__dirname, "tmp.txt"),
  //   JSON.stringify(data.dataY),
  //    // @ts-ignore
  //   function (err) {
  //     if (err) {
  //       return console.log(err);
  //     }
  //     console.log("The file was saved!");
  //   }
  // ); // Orfs.writeFileSync('/tmp/test-sync', 'Hey there!');
};

function testMemoryUsage() {
  // const arr = [1, 2, 3, 4, 5, 6, 9, 7, 8, 9, 10];
  // const arr = Array(1e7).fill(1e3);
  // arr.reverse();
  const used = process.memoryUsage().heapUsed / 1024 / 1024;
  console.log(`The script uses approximately ${Math.round(used * 100) / 100} MB`);
}


test1D();
test2D();
// testMemoryUsage();
