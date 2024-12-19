import { TeosBase } from '@greenroom-robotics/gsw_node';
import { describe, expect } from '@jest/globals';

describe('TeosBase', () => {
  const teosBase = new TeosBase();

  describe('gsw_z_from_p', () => {
    it('matches demo value 1', () => {
      const p = 10;
      const lat = 4;
      const expectedResult = -0.099445834469453 * 1.0e2;
      const result = teosBase.gsw_z_from_p(p, lat, 0, 0);

      expect(result).toBeCloseTo(expectedResult);
    });

    it('matches demo value 2', () => {
      const p = 50;
      const lat = 4;
      const expectedResult = -0.49718089701255 * 1.0e2;
      const result = teosBase.gsw_z_from_p(p, lat, 0.0, 0.0);

      expect(result).toBeCloseTo(expectedResult);
    });
  });

  describe('Calculate salinity', () => {
    it('matches demo value 1', () => {
      const C = 34.5487;
      const t = 28.7856;
      const p = 10;
      const expectedResult = 20.009869599086951;
      const result = teosBase.gsw_sp_from_c(C, t, p);

      expect(result).toBeCloseTo(expectedResult);
    });
  });
});
