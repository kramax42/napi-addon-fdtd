import fs from 'fs';
import path from 'path';

import addon from './index';



const test1D = () => {

  const condition = [1, 10, 1];
  const eps = [1, 1.2];
  const sigma = [0, 0.04];
  const srcPosition = [0.4, 0.8];

  let data = addon.getData2D(condition, true, eps, 2, srcPosition, sigma);

  for (let j = 0; j < 50; ++j) {
    data = addon.getData2D(condition, false, eps, 2, srcPosition, sigma);
  }

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
  const eps = Array(epsSize*2).fill(0).map(_ => Math.random()*10)

  const condition = [1, 10]

  let reload = true;
  let data = addon.getData3D(condition, reload, [4,5,6,7], 2, 0);

  reload = false;
  for (let j = 0; j < 150; ++j) {
    //eps, epsSize
    data = addon.getData3D(condition, reload, [4,5,6,7], 2, 0);
  }

  
  fs.writeFileSync(
    path.resolve(__dirname, "tmp.txt"),
    JSON.stringify(data.dataY),
     // @ts-ignore
    function (err) {
      if (err) {
        return console.log(err);
      }
      console.log("The file was saved!");
    }
  ); // Orfs.writeFileSync('/tmp/test-sync', 'Hey there!');
};

function testMemoryUsage() {
  // const arr = [1, 2, 3, 4, 5, 6, 9, 7, 8, 9, 10];
  // const arr = Array(1e7).fill(1e3);
  // arr.reverse();
  const used = process.memoryUsage().heapUsed / 1024 / 1024;
  console.log(`The script uses approximately ${Math.round(used * 100) / 100} MB`);
}


test1D();
// test2D();
// testMemoryUsage();