import addon from '@greenroom-robotics/gsw_cpp/build/Release/gsw-node.node'
import {TeosSea as TeosSeaClass} from './gsw-node-classes'

export const TeosSea: typeof TeosSeaClass = addon.TeosSea