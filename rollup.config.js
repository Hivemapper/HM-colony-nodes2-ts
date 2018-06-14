import typescript from 'rollup-plugin-typescript2';
import pkg from './package.json';

export default {
  input: 'src/export.ts',
  output: {
    file: 'dist/export.cjs.js',
    format: 'cjs',
    sourcemap: true,
  },
  plugins: [ typescript() ],
  external: Object.keys(pkg.dependencies),
};
