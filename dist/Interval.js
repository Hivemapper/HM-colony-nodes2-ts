import { S2 } from "./S2";
export class Interval {
    constructor(lo, hi) {
        this.lo = S2.toDecimal(lo);
        this.hi = S2.toDecimal(hi);
    }
    toString() {
        return "[" + this.lo.toString() + ", " + this.hi.toString() + "]";
    }
    /**
     * Return true if two intervals contains the same set of points.
     */
    equals(that) {
        if (typeof (that) === typeof (this)) {
            return this.lo.eq(that.lo) && this.hi.eq(that.hi);
        }
        return false;
    }
}
//# sourceMappingURL=Interval.js.map