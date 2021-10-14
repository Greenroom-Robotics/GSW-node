import * as nbind from '@mcesystems/nbind';
import * as LibTypes from './lib-types';

const lib = nbind.init<typeof LibTypes>(__dirname).lib;

export const TeosSea = lib.TeosSea;

export const TeosBase = lib.TeosBase;

export const TeosIce = lib.TeosIce;


