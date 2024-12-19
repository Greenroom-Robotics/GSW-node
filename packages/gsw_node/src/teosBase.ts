import addon from '@greenroom-robotics/gsw_cpp/build/Release/gsw-node.node';
import { TeosBase as TeosBaseClass } from './gsw-node-classes';

export const TeosBase: typeof TeosBaseClass = addon.TeosBase;
