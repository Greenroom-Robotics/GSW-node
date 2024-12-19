import { TeosIce } from '../src';
import { describe, expect } from '@jest/globals';

describe('TeosIce', () => {
  const teosIce = new TeosIce();

  describe('gsw_cp_ice', () => {
    it('matches demo value 1', () => {
      const t = -10.7856;
      const p = 10;
      const expectedResult = 2.017314262094657 * 1.0e3;
      const result = teosIce.gsw_cp_ice(t, p);

      expect(result).toBeCloseTo(expectedResult);
    });

    it('matches demo value 2', () => {
      const t = -13.4329;
      const p = 50;
      const expectedResult = 1.997830122682709 * 1.0e3;
      const result = teosIce.gsw_cp_ice(t, p);

      expect(result).toBeCloseTo(expectedResult);
    });
  });
});
