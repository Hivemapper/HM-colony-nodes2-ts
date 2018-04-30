import { Interval } from "./Interval";
import { S2 } from "./S2";
import * as decimal from 'decimal.js';
/**
 * An R1Interval represents a closed interval on a unit circle (also known as a
 * 1-dimensional sphere). It is capable of representing the empty interval
 * (containing no points), the full interval (containing all points), and
 * zero-length intervals (containing a single point).
 *
 *  Points are represented by the angle they make with the positive x-axis in
 * the range [-Pi, Pi]. An interval is represented by its lower and upper bounds
 * (both inclusive, since the interval is closed). The lower bound may be
 * greater than the upper bound, in which case the interval is "inverted" (i.e.
 * it passes through the point (-1, 0)).
 *
 *  Note that the point (-1, 0) has two valid representations, Pi and -Pi. The
 * normalized representation of this point internally is Pi, so that endpoints
 * of normal intervals are in the range (-Pi, Pi]. However, we take advantage of
 * the point -Pi to construct two special intervals: the Full() interval is
 * [-Pi, Pi], and the Empty() interval is [Pi, -Pi].
 *
 */
export class R1Interval extends Interval {
    /** Return true if the interval is empty, i.e. it contains no points. */
    isEmpty() {
        return this.lo.gt(this.hi);
    }
    getCenter() {
        return this.lo.plus(this.hi).dividedBy(2);
    }
    getLength() {
        return this.hi.minus(this.lo);
    }
    contains(_p) {
        const p = S2.toDecimal(_p);
        return p.gte(this.lo) && p.lte(this.hi);
    }
    /** Return true if the interior of the interval contains the point 'p'. */
    interiorContains(_p) {
        const p = S2.toDecimal(_p);
        return p.gt(this.lo) && p.lt(this.hi);
    }
    /**
     * Return true if the interval contains the given interval 'y'. Works for
     * empty, full, and singleton intervals.
     */
    containsI(y) {
        if (y.isEmpty()) {
            return true;
        }
        return y.lo.gte(this.lo) && y.hi.lte(this.hi);
    }
    interiorContainsI(y) {
        if (y.isEmpty()) {
            return true;
        }
        return y.lo.gt(this.lo) && y.hi.lt(this.hi);
    }
    /**
     * Return true if this interval intersects the given interval, i.e. if they
     * have any points in common.
     */
    intersects(y) {
        if (this.lo.lte(y.lo)) {
            return y.lo.lte(this.hi) && y.lo.lte(y.hi);
        }
        else {
            return this.lo.lte(y.hi) && this.lo.lte(this.hi);
        }
    }
    /**
     * Return true if the interior of this interval intersects any point of the
     * given interval (including its boundary).
     */
    interiorIntersects(y) {
        return y.lo.lt(this.hi) && this.lo.lt(y.hi) && this.lo.lt(this.hi) && y.lo.lte(y.hi);
    }
    /** Expand the interval so that it contains the given point "p". */
    addPoint(_p) {
        const p = S2.toDecimal(_p);
        if (this.isEmpty()) {
            return R1Interval.fromPoint(p);
        }
        else if (p.lt(this.lo)) {
            return new R1Interval(p, this.hi);
        }
        else if (p.gt(this.hi)) {
            return new R1Interval(this.lo, p);
        }
        else {
            return new R1Interval(this.lo, this.hi);
        }
    }
    /**
     * Return an interval that contains all points with a distance "radius" of a
     * point in this interval. Note that the expansion of an empty interval is
     * always empty.
     */
    expanded(_radius) {
        const radius = S2.toDecimal(_radius);
        // assert (radius >= 0);
        if (this.isEmpty()) {
            return this;
        }
        return new R1Interval(this.lo.minus(radius), this.hi.plus(radius));
    }
    /**
     * Return the smallest interval that contains this interval and the given
     * interval "y".
     */
    union(y) {
        if (this.isEmpty()) {
            return y;
        }
        if (y.isEmpty()) {
            return this;
        }
        return new R1Interval(decimal.Decimal.min(this.lo, y.lo), decimal.Decimal.max(this.hi, y.hi));
    }
    /**
     * Return the intersection of this interval with the given interval. Empty
     * intervals do not need to be special-cased.
     */
    intersection(y) {
        return new R1Interval(decimal.Decimal.max(this.lo, y.lo), decimal.Decimal.min(this.hi, y.hi));
    }
    /**
     * Return true if the length of the symmetric difference between the two
     * intervals is at most the given tolerance.
     */
    approxEquals(y, maxError = 1e-15) {
        if (this.isEmpty()) {
            return y.getLength().lte(maxError);
        }
        if (y.isEmpty()) {
            return this.getLength().lte(maxError);
        }
        return y.lo.minus(this.lo).abs()
            .plus(y.hi.minus(this.hi).abs())
            .lte(maxError);
    }
    static empty() {
        return new R1Interval(1, 0);
    }
    static fromPoint(p) {
        return new R1Interval(p, p);
    }
    /**
     * Convenience method to construct the minimal interval containing the two
     * given points. This is equivalent to starting with an empty interval and
     * calling AddPoint() twice, but it is more efficient.
     */
    static fromPointPair(_p1, _p2) {
        const p1 = S2.toDecimal(_p1);
        const p2 = S2.toDecimal(_p2);
        if (p1.lte(p2)) {
            return new R1Interval(p1, p2);
        }
        else {
            return new R1Interval(p2, p1);
        }
    }
}
//# sourceMappingURL=R1Interval.js.map