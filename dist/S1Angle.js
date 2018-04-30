import * as decimal from 'decimal.js';
import { S2 } from "./S2";
export class S1Angle {
    constructor(radians) {
        this.radians = new decimal.Decimal(radians);
    }
    degrees() {
        return S2.toDecimal(this.radians).times((180 / Math.PI));
    }
    //
    // public long e5() {
    //   return Math.round(degrees() * 1e5);
    // }
    //
    // public long e6() {
    //   return Math.round(degrees() * 1e6);
    // }
    //
    // public long e7() {
    //   return Math.round(degrees() * 1e7);
    // }
    /**
     * Return the angle between two points, which is also equal to the distance
     * between these points on the unit sphere. The points do not need to be
     * normalized.
     */
    static fromPoints(x, y) {
        return new S1Angle(x.angle(y));
    }
    lessThan(that) {
        return this.radians.lt(that.radians);
    }
    greaterThan(that) {
        return this.radians.gt(that.radians);
    }
    lessOrEquals(that) {
        return this.radians.lte(that.radians);
    }
    greaterOrEquals(that) {
        return this.radians.gte(that.radians);
    }
    static max(left, right) {
        return right.greaterThan(left) ? right : left;
    }
    static min(left, right) {
        return right.greaterThan(left) ? left : right;
    }
    static degrees(degrees) {
        let d = new decimal.Decimal(degrees);
        return new S1Angle(d.times(Math.PI / 180));
    }
    //
    // public static S1Angle e5(long e5) {
    //   return degrees(e5 * 1e-5);
    // }
    //
    // public static S1Angle e6(long e6) {
    //   // Multiplying by 1e-6 isn't quite as accurate as dividing by 1e6,
    //   // but it's about 10 times faster and more than accurate enough.
    //   return degrees(e6 * 1e-6);
    // }
    //
    // public static S1Angle e7(long e7) {
    //   return degrees(e7 * 1e-7);
    // }
    /**
     * Writes the angle in degrees with a "d" suffix, e.g. "17.3745d". By default
     * 6 digits are printed; this can be changed using setprecision(). Up to 17
     * digits are required to distinguish one angle from another.
     */
    toString() {
        return this.degrees() + "d";
    }
    compareTo(that) {
        return this.radians < that.radians ? -1 : this.radians > that.radians ? 1 : 0;
    }
}
//# sourceMappingURL=S1Angle.js.map