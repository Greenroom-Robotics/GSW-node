import { TeosSea } from '../src';
import { describe, expect } from '@jest/globals';

describe('TeosSea', () => {
  const teosSea = new TeosSea();

  describe('gsw_C_from_SP', () => {
    it('matches demo value 1', () => {
      const SP = 34.5487;
      const t = 28.7856;
      const p = 10;
      const expectedResult = 56.412599581571186;
      const result = teosSea.gsw_c_from_sp(SP, t, p);

      expect(result).toBeCloseTo(expectedResult);
    });

    it('matches demo value 2', () => {
      const SP = 34.7275;
      const t = 28.4329;
      const p = 50;
      const expectedResult = 56.316185602699953;
      const result = teosSea.gsw_c_from_sp(SP, t, p);

      expect(result).toBeCloseTo(expectedResult);
    });
  });
});
