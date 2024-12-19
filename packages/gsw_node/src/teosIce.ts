import addon from '@greenroom-robotics/gsw_cpp/build/Release/gsw-node.node';
import { TeosIce as TeosIceClass } from './gsw-node-classes';

export const TeosIce: typeof TeosIceClass = addon.TeosIce;
