import { S2 } from "./S2";
var Interval = /** @class */ (function () {
    function Interval(lo, hi) {
        this.lo = S2.toDecimal(lo);
        this.hi = S2.toDecimal(hi);
    }
    Interval.prototype.toString = function () {
        return "[" + this.lo.toString() + ", " + this.hi.toString() + "]";
    };
    /**
     * Return true if two intervals contains the same set of points.
     */
    Interval.prototype.equals = function (that) {
        if (typeof (that) === typeof (this)) {
            return this.lo.eq(that.lo) && this.hi.eq(that.hi);
        }
        return false;
    };
    return Interval;
}());
export { Interval };
//# sourceMappingURL=Interval.js.map