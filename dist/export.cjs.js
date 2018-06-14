'use strict';

Object.defineProperty(exports, '__esModule', { value: true });

function _interopDefault (ex) { return (ex && (typeof ex === 'object') && 'default' in ex) ? ex['default'] : ex; }

var decimal = require('decimal.js');
var Long = _interopDefault(require('long'));

/**
 * R2Vector represents a vector in the two-dimensional space. It defines the
 * basic geometrical operations for 2D vectors, e.g. cross product, addition,
 * norm, comparison etc.
 *
 */
var R2Vector = /** @class */ (function () {
    function R2Vector(_x, _y) {
        this._x = new decimal.Decimal(_x);
        this._y = new decimal.Decimal(_y);
        // this._x = new Decimal(_x) as decimal.Decimal;
        // this._y = new Decimal(_y) as decimal.Decimal;
    }
    Object.defineProperty(R2Vector.prototype, "x", {
        get: function () {
            return this._x;
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(R2Vector.prototype, "y", {
        get: function () {
            return this._y;
        },
        enumerable: true,
        configurable: true
    });
    R2Vector.prototype.get = function (index) {
        if (index > 1) {
            throw new Error("Index out fo bounds error " + index);
        }
        return index == 0 ? this._x : this._y;
    };
    R2Vector.fromPointFace = function (p, face) {
        return p.toR2Vector(face);
    };
    R2Vector.add = function (p1, p2) {
        return new R2Vector(p1._x.plus(p2._x), p1._y.plus(p2._y));
    };
    R2Vector.mul = function (p, _m) {
        var m = new decimal.Decimal(_m);
        return new R2Vector(m.times(p._x), m.times(p._y));
    };
    R2Vector.prototype.norm2 = function () {
        return this.x.pow(2).plus(this.y.pow(2));
    };
    R2Vector.dotProd = function (p1, p2) {
        return p1.x.times(p2.x).plus(p1.y.times(p2.y));
    };
    R2Vector.prototype.dotProd = function (that) {
        return R2Vector.dotProd(this, that);
    };
    R2Vector.prototype.crossProd = function (that) {
        return this.x.times(that.y).minus(this.y.times(that.x));
    };
    R2Vector.prototype.lessThan = function (vb) {
        if (this.x.lt(vb.x)) {
            return true;
        }
        if (vb.x.lt(this.x)) {
            return false;
        }
        if (this.y.lt(vb.y)) {
            return true;
        }
        return false;
    };
    //
    // @Override
    // public boolean equals(Object that) {
    //   if (!(that instanceof R2Vector)) {
    //     return false;
    //   }
    //   R2Vector thatPoint = (R2Vector) that;
    //   return this.x == thatPoint.x && this.y == thatPoint.y;
    // }
    // /**
    //  * Calcualates hashcode based on stored coordinates. Since we want +0.0 and
    //  * -0.0 to be treated the same, we ignore the sign of the coordinates.
    //  */
    // @Override
    // public int hashCode() {
    //   long value = 17;
    //   value += 37 * value + Double.doubleToLongBits(Math.abs(x));
    //   value += 37 * value + Double.doubleToLongBits(Math.abs(y));
    //   return (int) (value ^ (value >>> 32));
    // }
    //
    R2Vector.fromSTVector = function (stVector) {
        return new R2Vector(R2Vector.singleStTOUV(stVector.x), R2Vector.singleStTOUV(stVector.y));
    };
    // from S2Projections.stToUV (QUADRATIC)
    R2Vector.singleStTOUV = function (_s) {
        var s = S2.toDecimal(_s);
        if (s.gte(0)) {
            return S2.toDecimal(1)
                .dividedBy(3)
                .times(s.plus(1).pow(2).minus(1));
            // return (1 / 3.) * ((1 + s) * (1 + s) - 1);
        }
        else {
            return S2.toDecimal(1)
                .dividedBy(3)
                .times(S2.toDecimal(1)
                .minus(S2.toDecimal(1).minus(s).pow(2)));
            // return (1 / 3.) * (1 - (1 - s) * (1 - s));
        }
    };
    R2Vector.singleUVToST = function (_x) {
        var x = S2.toDecimal(_x);
        if (x.gte(0)) {
            return decimal.Decimal.sqrt(x.times(3).plus(1)).minus(1);
        }
        else {
            return S2.toDecimal(1)
                .minus(decimal.Decimal.sqrt(S2.toDecimal(1).minus(x.times(3))));
        }
    };
    /**
     * To be used only if this vector is representing uv.
     * @param face
     * @returns {S2Point}
     */
    R2Vector.prototype.toPoint = function (face) {
        switch (face) {
            case 0:
                return new S2Point(1, this.x, this.y);
            case 1:
                return new S2Point(this.x.neg(), 1, this.y);
            case 2:
                return new S2Point(this.x.neg(), this.y.neg(), 1);
            case 3:
                return new S2Point(-1, this.y.neg(), this.x.neg());
            case 4:
                return new S2Point(this.y, -1, this.x.neg());
            default:
                return new S2Point(this.y, this.x, -1);
        }
    };
    R2Vector.prototype.toSt = function (which) {
        return which == 0 ? R2Vector.singleUVToST(this.x) : R2Vector.singleUVToST(this.y);
    };
    R2Vector.prototype.toString = function () {
        return "(" + this.x.toString() + ", " + this.y.toString() + ")";
    };
    return R2Vector;
}());

/*
 * Copyright 2006 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
///re
/**
 * An S2Point represents a point on the unit sphere as a 3D vector. Usually
 * points are normalized to be unit length, but some methods do not require
 * this.
 *
 */
var S2Point = /** @class */ (function () {
    function S2Point(x, y, z) {
        this.x = new decimal.Decimal(x);
        this.y = new decimal.Decimal(y);
        this.z = new decimal.Decimal(z);
        // this.y = typeof(y) === 'number'?new decimal.Decimal(y):y as decimal.Decimal;
        // this.z = typeof(z) === 'number'?new decimal.Decimal(z):z as decimal.Decimal;
    }
    S2Point.minus = function (p1, p2) {
        return S2Point.sub(p1, p2);
    };
    S2Point.neg = function (p) {
        return new S2Point(p.x.negated(), p.y.negated(), p.z.negated());
    };
    S2Point.prototype.norm2 = function () {
        return this.x.pow(2).plus(this.y.pow(2)).plus(this.z.pow(2));
    };
    S2Point.prototype.norm = function () {
        return this.norm2().sqrt();
    };
    S2Point.crossProd = function (p1, p2) {
        return new S2Point(p1.y.times(p2.z).minus(p1.z.times(p2.y)), p1.z.times(p2.x).minus(p1.x.times(p2.z)), 
        // p1.z * p2.x - p1.x * p2.z,
        p1.x.times(p2.y).minus(p1.y.times(p2.x))
        // p1.x * p2.y - p1.y * p2.x
        );
    };
    S2Point.add = function (p1, p2) {
        return new S2Point(p1.x.add(p2.x), p1.y.add(p2.y), p1.z.add(p2.z));
    };
    S2Point.sub = function (p1, p2) {
        return new S2Point(p1.x.sub(p2.x), p1.y.sub(p2.y), p1.z.sub(p2.z));
    };
    S2Point.prototype.dotProd = function (that) {
        return this.x.times(that.x).plus(this.y.times(that.y)).plus(this.z.times(that.z));
    };
    S2Point.mul = function (p, m) {
        var mD = new decimal.Decimal(m);
        return new S2Point(mD.times(p.x), mD.times(p.y), mD.times(p.z));
    };
    S2Point.div = function (p, m) {
        return new S2Point(p.x.div(m), p.y.div(m), p.z.div(m));
    };
    /** return a vector orthogonal to this one */
    S2Point.prototype.ortho = function () {
        var k = this.largestAbsComponent();
        var temp;
        if (k == 1) {
            temp = new S2Point(1, 0, 0);
        }
        else if (k == 2) {
            temp = new S2Point(0, 1, 0);
        }
        else {
            temp = new S2Point(0, 0, 1);
        }
        return S2Point.normalize(S2Point.crossProd(this, temp));
    };
    /** Return the index of the largest component fabs */
    S2Point.prototype.largestAbsComponent = function () {
        var temp = S2Point.fabs(this);
        if (temp.x.greaterThan(temp.y)) {
            if (temp.x.greaterThan(temp.z)) {
                return 0;
            }
            else {
                return 2;
            }
        }
        else {
            if (temp.y.greaterThan(temp.z)) {
                return 1;
            }
            else {
                return 2;
            }
        }
    };
    S2Point.fabs = function (p) {
        return new S2Point(p.x.abs(), p.y.abs(), p.z.abs());
    };
    S2Point.normalize = function (p) {
        var norm = p.norm();
        if (!norm.eq(0)) {
            norm = S2.toDecimal(1).dividedBy(norm);
        }
        return S2Point.mul(p, norm);
    };
    S2Point.prototype.axis = function (axis) {
        return (axis == 0) ? this.x : (axis == 1) ? this.y : this.z;
    };
    /** Return the angle between two vectors in radians */
    S2Point.prototype.angle = function (va) {
        return decimal.Decimal.atan2(S2Point.crossProd(this, va).norm(), this.dotProd(va));
    };
    /**
     * Compare two vectors, return true if all their components are within a
     * difference of margin.
     */
    S2Point.prototype.aequal = function (that, margin) {
        return this.x.minus(that.x).abs().lessThan(margin) &&
            this.y.minus(that.y).abs().lessThan(margin) &&
            this.z.minus(that.z).abs().lessThan(margin);
    };
    S2Point.prototype.equals = function (that) {
        if (!(that instanceof S2Point)) {
            return false;
        }
        return this.x.eq(that.x) && this.y.eq(that.y) && this.z.eq(that.z);
    };
    S2Point.prototype.lessThan = function (vb) {
        if (this.x.lt(vb.x)) {
            return true;
        }
        if (vb.x.lt(this.x)) {
            return false;
        }
        if (this.y.lt(vb.y)) {
            return true;
        }
        if (vb.y.lt(this.y)) {
            return false;
        }
        if (this.z.lt(vb.z)) {
            return true;
        }
        return false;
    };
    S2Point.prototype.compareTo = function (other) {
        return (this.lessThan(other) ? -1 : (this.equals(other) ? 0 : 1));
    };
    S2Point.prototype.toFace = function () {
        var face = this.largestAbsComponent();
        if (this.axis(face).lt(0)) {
            face += 3;
        }
        return face;
    };
    S2Point.prototype.toR2Vector = function (face) {
        if (face === void 0) { face = this.toFace(); }
        var u;
        var v;
        switch (face) {
            case 0:
                u = this.y.div(this.x);
                v = this.z.div(this.x);
                break;
            case 1:
                u = this.x.neg().div(this.y);
                v = this.z.div(this.y);
                break;
            case 2:
                u = this.x.neg().div(this.z);
                v = this.y.neg().div(this.z);
                break;
            case 3:
                u = this.z.div(this.x);
                v = this.y.div(this.x);
                break;
            case 4:
                u = this.z.div(this.y);
                v = this.x.neg().div(this.y);
                break;
            case 5:
                u = this.y.neg().div(this.z);
                v = this.x.neg().div(this.z);
                break;
            default:
                throw new Error('Invalid face');
        }
        return new R2Vector(u, v);
    };
    S2Point.prototype.toString = function () {
        return "Point(" + this.x.toNumber() + ", " + this.y.toNumber() + ", " + this.z.toNumber() + ")";
    };
    return S2Point;
}());

/**
 * Defines an area or a length cell metric.
 */
var S2Metric$$1 = /** @class */ (function () {
    /**
     * Defines a cell metric of the given dimension (1 == length, 2 == area).
     */
    function S2Metric$$1(_dim, _deriv) {
        this._dim = S2.toDecimal(_dim).toNumber();
        this._deriv = S2.toDecimal(_deriv);
    }
    S2Metric$$1.prototype.deriv = function () {
        return this._deriv;
    };
    S2Metric$$1.prototype.dim = function () {
        return this._dim;
    };
    /** Return the value of a metric for cells at the given level. */
    S2Metric$$1.prototype.getValue = function (level) {
        var scaleFactor = this.dim() * (1 - level);
        return this.deriv().times(Math.pow(2, scaleFactor)).toNumber();
    };
    /**
     * Return the level at which the metric has approximately the given value.
     * For example, S2::kAvgEdge.GetClosestLevel(0.1) returns the level at which
     * the average cell edge length is approximately 0.1. The return value is
     * always a valid level.
     */
    S2Metric$$1.prototype.getClosestLevel = function (/*double*/ value) {
        return this.getMinLevel(S2.M_SQRT2 * value);
    };
    /**
     * Return the minimum level such that the metric is at most the given value,
     * or S2CellId::kMaxLevel if there is no such level. For example,
     * S2::kMaxDiag.GetMinLevel(0.1) returns the minimum level such that all
     * cell diagonal lengths are 0.1 or smaller. The return value is always a
     * valid level.
     */
    S2Metric$$1.prototype.getMinLevel = function (value /*double*/) {
        if (value <= 0) {
            return S2.MAX_LEVEL;
        }
        // This code is equivalent to computing a floating-point "level"
        // value and rounding up.
        var exponent = S2.exp(value / ((1 << this.dim()) * this.deriv().toNumber()));
        var level = Math.max(0, Math.min(S2.MAX_LEVEL, -((exponent - 1) >> (this.dim() - 1))));
        // assert (level == S2CellId.MAX_LEVEL || getValue(level) <= value);
        // assert (level == 0 || getValue(level - 1) > value);
        return level;
    };
    /**
     * Return the maximum level such that the metric is at least the given
     * value, or zero if there is no such level. For example,
     * S2.kMinWidth.GetMaxLevel(0.1) returns the maximum level such that all
     * cells have a minimum width of 0.1 or larger. The return value is always a
     * valid level.
     */
    S2Metric$$1.prototype.getMaxLevel = function (_value /*double*/) {
        var value = S2.toDecimal(_value).toNumber();
        if (value <= 0) {
            return S2.MAX_LEVEL;
        }
        // This code is equivalent to computing a floating-point "level"
        // value and rounding down.
        var exponent = S2.exp((1 << this.dim()) * this.deriv().toNumber() / value);
        var level = Math.max(0, Math.min(S2.MAX_LEVEL, ((exponent - 1) >> (this.dim() - 1))));
        // assert (level == 0 || getValue(level) >= value);
        // assert (level == S2CellId.MAX_LEVEL || getValue(level + 1) < value);
        return level;
    };
    return S2Metric$$1;
}());

var S2 = /** @class */ (function () {
    function S2() {
    }
    S2.IEEEremainder = function (_f1, _f2) {
        var f1 = S2.toDecimal(_f1);
        var f2 = S2.toDecimal(_f2);
        var r = f1.mod(f2);
        if (r.isNaN() || r.eq(f2) || r.lessThanOrEqualTo(f2.abs().dividedBy(2))) {
            return r;
        }
        else {
            return (f1.gte(0) ? S2.toDecimal(1) : S2.toDecimal(-1)).times(r.minus(f2));
        }
    };
    /**
     * Return true if the given point is approximately unit length (this is mainly
     * useful for assertions).
     */
    S2.isUnitLength = function (p) {
        return p.norm2().minus(1).abs().lte(1e-15);
    };
    /**
     * If v is non-zero, return an integer {@code exp} such that
     * {@code (0.5 <= |v|*2^(-exp) < 1)}. If v is zero, return 0.
     *
     * <p>Note that this arguably a bad definition of exponent because it makes
     * {@code exp(9) == 4}. In decimal this would be like saying that the
     * exponent of 1234 is 4, when in scientific 'exponent' notation 1234 is
     * {@code 1.234 x 10^3}.
     *
     * TODO(dbeaumont): Replace this with "DoubleUtils.getExponent(v) - 1" ?
     */
    S2.exp = function (v /*double*/) {
        if (v == 0) {
            return 0;
        }
        // IT should always be ((int)log(2,v))+1;
        var start = Math.floor(Math.log(v) / Math.log(2));
        for (var i = start; i < start + 10; i++) {
            var curVal = Math.abs(v) * Math.pow(2, -i);
            if (curVal >= 0.5 && curVal < 1) {
                return i;
            }
        }
        throw new Error('method not written yet');
        // return (int)((S2.EXPONENT_MASK & bits) >> S2.EXPONENT_SHIFT) - 1022;
    };
    /**
     * Return a vector "c" that is orthogonal to the given unit-length vectors "a"
     * and "b". This function is similar to a.CrossProd(b) except that it does a
     * better job of ensuring orthogonality when "a" is nearly parallel to "b",
     * and it returns a non-zero result even when a == b or a == -b.
     *
     *  It satisfies the following properties (RCP == RobustCrossProd):
     *
     *  (1) RCP(a,b) != 0 for all a, b (2) RCP(b,a) == -RCP(a,b) unless a == b or
     * a == -b (3) RCP(-a,b) == -RCP(a,b) unless a == b or a == -b (4) RCP(a,-b)
     * == -RCP(a,b) unless a == b or a == -b
     */
    S2.robustCrossProd = function (a, b) {
        // The direction of a.CrossProd(b) becomes unstable as (a + b) or (a - b)
        // approaches zero. This leads to situations where a.CrossProd(b) is not
        // very orthogonal to "a" and/or "b". We could fix this using Gram-Schmidt,
        // but we also want b.RobustCrossProd(a) == -b.RobustCrossProd(a).
        //
        // The easiest fix is to just compute the cross product of (b+a) and (b-a).
        // Given that "a" and "b" are unit-length, this has good orthogonality to
        // "a" and "b" even if they differ only in the lowest bit of one component.
        // assert (isUnitLength(a) && isUnitLength(b));
        var x = S2Point.crossProd(S2Point.add(b, a), S2Point.sub(b, a));
        if (!x.equals(new S2Point(0, 0, 0))) {
            return x;
        }
        // The only result that makes sense mathematically is to return zero, but
        // we find it more convenient to return an arbitrary orthogonal vector.
        return a.ortho();
    };
    /**
     * Return the area of triangle ABC. The method used is about twice as
     * expensive as Girard's formula, but it is numerically stable for both large
     * and very small triangles. The points do not need to be normalized. The area
     * is always positive.
     *
     *  The triangle area is undefined if it contains two antipodal points, and
     * becomes numerically unstable as the length of any edge approaches 180
     * degrees.
     */
    S2.area = function (a, b, c) {
        // This method is based on l'Huilier's theorem,
        //
        // tan(E/4) = sqrt(tan(s/2) tan((s-a)/2) tan((s-b)/2) tan((s-c)/2))
        //
        // where E is the spherical excess of the triangle (i.e. its area),
        // a, b, c, are the side lengths, and
        // s is the semiperimeter (a + b + c) / 2 .
        //
        // The only significant source of error using l'Huilier's method is the
        // cancellation error of the terms (s-a), (s-b), (s-c). This leads to a
        // *relative* error of about 1e-16 * s / min(s-a, s-b, s-c). This compares
        // to a relative error of about 1e-15 / E using Girard's formula, where E is
        // the true area of the triangle. Girard's formula can be even worse than
        // this for very small triangles, e.g. a triangle with a true area of 1e-30
        // might evaluate to 1e-5.
        //
        // So, we prefer l'Huilier's formula unless dmin < s * (0.1 * E), where
        // dmin = min(s-a, s-b, s-c). This basically includes all triangles
        // except for extremely long and skinny ones.
        //
        // Since we don't know E, we would like a conservative upper bound on
        // the triangle area in terms of s and dmin. It's possible to show that
        // E <= k1 * s * sqrt(s * dmin), where k1 = 2*sqrt(3)/Pi (about 1).
        // Using this, it's easy to show that we should always use l'Huilier's
        // method if dmin >= k2 * s^5, where k2 is about 1e-2. Furthermore,
        // if dmin < k2 * s^5, the triangle area is at most k3 * s^4, where
        // k3 is about 0.1. Since the best case error using Girard's formula
        // is about 1e-15, this means that we shouldn't even consider it unless
        // s >= 3e-4 or so.
        // We use volatile doubles to force the compiler to truncate all of these
        // quantities to 64 bits. Otherwise it may compute a value of dmin > 0
        // simply because it chose to spill one of the intermediate values to
        // memory but not one of the others.
        var sa = b.angle(c);
        var sb = c.angle(a);
        var sc = a.angle(b);
        var s = sa.plus(sb).plus(sc).times(0.5);
        // 0.5 * (sa + sb + sc);
        if (s.gte(3e-4)) {
            // Consider whether Girard's formula might be more accurate.
            var s2 = s.pow(2);
            var dmin = s.minus(decimal.Decimal.max(sa, sb, sc));
            if (dmin.lt(s2.pow(2).times(s).times(1e-2))) {
                // This triangle is skinny enough to consider Girard's formula.
                var area = S2.girardArea(a, b, c);
                if (dmin.lt(s.times(area.times(0.1)))) {
                    return area;
                }
            }
        }
        // Use l'Huilier's formula.
        return S2.toDecimal(4)
            .times(decimal.Decimal.atan(decimal.Decimal.sqrt(decimal.Decimal.max(0.0, decimal.Decimal.tan(s.times(0.5))
            .times(decimal.Decimal.tan(s.minus(sa).times(0.5)))
            .times(decimal.Decimal.tan(s.minus(sb).times(0.5)))
            .times(decimal.Decimal.tan(s.minus(sc).times(0.5)))))));
    };
    /**
     * Return the area of the triangle computed using Girard's formula. This is
     * slightly faster than the Area() method above is not accurate for very small
     * triangles.
     */
    S2.girardArea = function (a, b, c) {
        // This is equivalent to the usual Girard's formula but is slightly
        // more accurate, faster to compute, and handles a == b == c without
        // a special case.
        var ab = S2Point.crossProd(a, b);
        var bc = S2Point.crossProd(b, c);
        var ac = S2Point.crossProd(a, c);
        return decimal.Decimal.max(0, ab.angle(ac)
            .minus(ab.angle(bc))
            .plus(bc.angle(ac)));
    };
    S2.toDecimal = function (value) {
        if (typeof (value) === 'number' || typeof (value) === 'string') {
            return new decimal.Decimal(value);
        }
        return value;
    };
    /**
     * Return true if the points A, B, C are strictly counterclockwise. Return
     * false if the points are clockwise or colinear (i.e. if they are all
     * contained on some great circle).
     *
     *  Due to numerical errors, situations may arise that are mathematically
     * impossible, e.g. ABC may be considered strictly CCW while BCA is not.
     * However, the implementation guarantees the following:
     *
     *  If SimpleCCW(a,b,c), then !SimpleCCW(c,b,a) for all a,b,c.
     *
     * In other words, ABC and CBA are guaranteed not to be both CCW
     */
    S2.simpleCCW = function (a, b, c) {
        // We compute the signed volume of the parallelepiped ABC. The usual
        // formula for this is (AxB).C, but we compute it here using (CxA).B
        // in order to ensure that ABC and CBA are not both CCW. This follows
        // from the following identities (which are true numerically, not just
        // mathematically):
        //
        // (1) x.CrossProd(y) == -(y.CrossProd(x))
        // (2) (-x).DotProd(y) == -(x.DotProd(y))
        return S2Point.crossProd(c, a).dotProd(b).gt(0);
    };
    /**
     *
     * Return true if edge AB crosses CD at a point that is interior to both
     * edges. Properties:
     *
     *  (1) SimpleCrossing(b,a,c,d) == SimpleCrossing(a,b,c,d) (2)
     * SimpleCrossing(c,d,a,b) == SimpleCrossing(a,b,c,d)
     */
    S2.simpleCrossing = function (a, b, c, d) {
        // We compute SimpleCCW() for triangles ACB, CBD, BDA, and DAC. All
        // of these triangles need to have the same orientation (CW or CCW)
        // for an intersection to exist. Note that this is slightly more
        // restrictive than the corresponding definition for planar edges,
        // since we need to exclude pairs of line segments that would
        // otherwise "intersect" by crossing two antipodal points.
        var ab = S2Point.crossProd(a, b);
        var cd = S2Point.crossProd(c, d);
        var acb = ab.dotProd(c).neg();
        var cbd = cd.dotProd(b).neg();
        var bda = ab.dotProd(d);
        var dac = cd.dotProd(a);
        return (acb.times(cbd).gt(0)) && (cbd.times(bda).gt(0)) && (bda.times(dac).gt(0));
    };
    S2.M_PI = Math.PI;
    S2.M_1_PI = 1.0 / Math.PI;
    S2.M_PI_2 = Math.PI / 2.0;
    S2.M_PI_4 = Math.PI / 4.0;
    S2.M_SQRT2 = Math.sqrt(2);
    S2.M_E = Math.E;
    // the axis directions are reversed).
    S2.SWAP_MASK = 0x01;
    S2.INVERT_MASK = 0x02;
    // Number of bits in the mantissa of a double.
    S2.EXPONENT_SHIFT = 52;
    // Mask to extract the exponent from a double.
    S2.EXPONENT_MASK = Long.fromString('0x7ff0000000000000', true, 16);
    /** Mapping from cell orientation + Hilbert traversal to IJ-index. */
    S2.POS_TO_ORIENTATION = [S2.SWAP_MASK, 0, 0, S2.INVERT_MASK + S2.SWAP_MASK];
    S2.POS_TO_IJ = [
        // 0 1 2 3
        [0, 1, 3, 2],
        [0, 2, 3, 1],
        [3, 2, 0, 1],
        [3, 1, 0, 2],
    ];
    S2.MAX_LEVEL = 30;
    S2.Metric = S2Metric$$1;
    return S2;
}());

var S1Angle = /** @class */ (function () {
    function S1Angle(radians) {
        this.radians = new decimal.Decimal(radians);
    }
    S1Angle.prototype.degrees = function () {
        return S2.toDecimal(this.radians).times((180 / Math.PI));
    };
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
    S1Angle.fromPoints = function (x, y) {
        return new S1Angle(x.angle(y));
    };
    S1Angle.prototype.lessThan = function (that) {
        return this.radians.lt(that.radians);
    };
    S1Angle.prototype.greaterThan = function (that) {
        return this.radians.gt(that.radians);
    };
    S1Angle.prototype.lessOrEquals = function (that) {
        return this.radians.lte(that.radians);
    };
    S1Angle.prototype.greaterOrEquals = function (that) {
        return this.radians.gte(that.radians);
    };
    S1Angle.max = function (left, right) {
        return right.greaterThan(left) ? right : left;
    };
    S1Angle.min = function (left, right) {
        return right.greaterThan(left) ? left : right;
    };
    S1Angle.degrees = function (degrees) {
        var d = new decimal.Decimal(degrees);
        return new S1Angle(d.times(Math.PI / 180));
    };
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
    S1Angle.prototype.toString = function () {
        return this.degrees() + "d";
    };
    S1Angle.prototype.compareTo = function (that) {
        return this.radians < that.radians ? -1 : this.radians > that.radians ? 1 : 0;
    };
    return S1Angle;
}());

/*! *****************************************************************************
Copyright (c) Microsoft Corporation. All rights reserved.
Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at http://www.apache.org/licenses/LICENSE-2.0

THIS CODE IS PROVIDED ON AN *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
MERCHANTABLITY OR NON-INFRINGEMENT.

See the Apache Version 2.0 License for specific language governing permissions
and limitations under the License.
***************************************************************************** */
/* global Reflect, Promise */

var extendStatics = Object.setPrototypeOf ||
    ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
    function (d, b) { for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p]; };

function __extends(d, b) {
    extendStatics(d, b);
    function __() { this.constructor = d; }
    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
}

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

var S1Interval = /** @class */ (function (_super) {
    __extends(S1Interval, _super);
    function S1Interval(lo, hi, checked) {
        if (checked === void 0) { checked = false; }
        var _this = _super.call(this, lo, hi) || this;
        if (!checked) {
            if (_this.lo.eq(-S2.M_PI) && !_this.hi.eq(S2.M_PI)) {
                _this.lo = S2.toDecimal(S2.M_PI);
            }
            if (_this.hi.eq(-S2.M_PI) && !_this.lo.eq(S2.M_PI)) {
                _this.hi = S2.toDecimal(S2.M_PI);
            }
        }
        return _this;
    }
    /**
     * An interval is valid if neither bound exceeds Pi in absolute value, and the
     * value -Pi appears only in the Empty() and Full() intervals.
     */
    S1Interval.prototype.isValid = function () {
        return this.lo.abs().lte(S2.M_PI) && this.hi.abs().lte(S2.M_PI)
            && !(this.lo.eq(-S2.M_PI) && !this.hi.eq(S2.M_PI))
            && !(this.hi.eq(-S2.M_PI) && !this.lo.eq(S2.M_PI));
        // return (Math.abs(this.lo) <= S2.M_PI && Math.abs(this.hi) <= S2.M_PI
        // && !(this.lo == -S2.M_PI && this.hi != S2.M_PI) && !(this.hi == -S2.M_PI && this.lo != S2.M_PI));
    };
    /** Return true if the interval contains all points on the unit circle. */
    S1Interval.prototype.isFull = function () {
        // console.log(this.hi.minus(this.lo).eq(2 * S2.M_PI));
        return this.hi.minus(this.lo).eq(2 * S2.M_PI);
    };
    /** Return true if the interval is empty, i.e. it contains no points. */
    S1Interval.prototype.isEmpty = function () {
        return this.lo.minus(this.hi).eq(2 * S2.M_PI);
    };
    /* Return true if this.lo > this.hi. (This is true for empty intervals.) */
    S1Interval.prototype.isInverted = function () {
        return this.lo.gt(this.hi);
    };
    /**
     * Return the midpoint of the interval. For full and empty intervals, the
     * result is arbitrary.
     */
    S1Interval.prototype.getCenter = function () {
        var center = this.lo.plus(this.hi).dividedBy(2);
        // let center = 0.5 * (this.lo + this.hi);
        if (!this.isInverted()) {
            return center;
        }
        // Return the center in the range (-Pi, Pi].
        return (center.lte(0)) ? (center.plus(S2.M_PI)) : (center.minus(S2.M_PI));
    };
    /**
     * Return the length of the interval. The length of an empty interval is
     * negative.
     */
    S1Interval.prototype.getLength = function () {
        var length = this.hi.minus(this.lo);
        if (length.gte(0)) {
            return length;
        }
        length = length.plus(2 * S2.M_PI);
        // Empty intervals have a negative length.
        return (length.gt(0)) ? length : S2.toDecimal(-1);
    };
    /**
     * Return the complement of the interior of the interval. An interval and its
     * complement have the same boundary but do not share any interior values. The
     * complement operator is not a bijection, since the complement of a singleton
     * interval (containing a single value) is the same as the complement of an
     * empty interval.
     */
    S1Interval.prototype.complement = function () {
        if (this.lo.eq(this.hi)) {
            return S1Interval.full(); // Singleton.
        }
        return new S1Interval(this.hi, this.lo, true); // Handles
        // empty and
        // full.
    };
    /** Return true if the interval (which is closed) contains the point 'p'. */
    S1Interval.prototype.contains = function (_p) {
        var p = S2.toDecimal(_p);
        // Works for empty, full, and singleton intervals.
        // assert (Math.abs(p) <= S2.M_PI);
        if (p.eq(-S2.M_PI)) {
            p = S2.toDecimal(S2.M_PI);
        }
        return this.fastContains(p);
    };
    /**
     * Return true if the interval (which is closed) contains the point 'p'. Skips
     * the normalization of 'p' from -Pi to Pi.
     *
     */
    S1Interval.prototype.fastContains = function (_p) {
        var p = S2.toDecimal(_p);
        if (this.isInverted()) {
            return (p.gte(this.lo) || p.lte(this.hi)) && !this.isEmpty();
        }
        else {
            return p.gte(this.lo) && p.lte(this.hi);
        }
    };
    /** Return true if the interior of the interval contains the point 'p'. */
    S1Interval.prototype.interiorContains = function (_p) {
        // Works for empty, full, and singleton intervals.
        // assert (Math.abs(p) <= S2.M_PI);
        var p = S2.toDecimal(_p);
        if (p.eq(-S2.M_PI)) {
            p = S2.toDecimal(S2.M_PI);
        }
        if (this.isInverted()) {
            return p.gt(this.lo) || p.lt(this.hi);
        }
        else {
            return (p.gt(this.lo) && p.lt(this.hi)) || this.isFull();
        }
    };
    /**
     * Return true if the interval contains the given interval 'y'. Works for
     * empty, full, and singleton intervals.
     */
    S1Interval.prototype.containsI = function (y) {
        // It might be helpful to compare the structure of these tests to
        // the simpler Contains(number) method above.
        if (this.isInverted()) {
            if (y.isInverted()) {
                return y.lo.gte(this.lo) && y.hi.lte(this.hi);
            }
            return (y.lo.gte(this.lo) || y.hi.lte(this.hi)) && !this.isEmpty();
        }
        else {
            if (y.isInverted()) {
                return this.isFull() || y.isEmpty();
            }
            return y.lo.gte(this.lo) && y.hi.lte(this.hi);
        }
    };
    /**
     * Returns true if the interior of this interval contains the entire interval
     * 'y'. Note that x.InteriorContains(x) is true only when x is the empty or
     * full interval, and x.InteriorContains(S1Interval(p,p)) is equivalent to
     * x.InteriorContains(p).
     */
    S1Interval.prototype.interiorContainsI = function (y) {
        if (this.isInverted()) {
            if (!y.isInverted()) {
                return this.lo.gt(this.lo) || y.hi.lt(this.hi);
            }
            return (y.lo.gt(this.lo) && y.hi.lt(this.hi)) || y.isEmpty();
        }
        else {
            if (y.isInverted()) {
                return this.isFull() || y.isEmpty();
            }
            return (y.lo.gt(this.lo) && y.hi.lt(this.hi)) || this.isFull();
        }
    };
    /**
     * Return true if the two intervals contain any points in common. Note that
     * the point +/-Pi has two representations, so the intervals [-Pi,-3] and
     * [2,Pi] intersect, for example.
     */
    S1Interval.prototype.intersects = function (y) {
        if (this.isEmpty() || y.isEmpty()) {
            return false;
        }
        if (this.isInverted()) {
            // Every non-empty inverted interval contains Pi.
            return y.isInverted() || y.lo.lte(this.hi) || y.hi.gte(this.lo);
        }
        else {
            if (y.isInverted()) {
                return y.lo.lte(this.hi) || y.hi.gte(this.lo);
            }
            return y.lo.lte(this.hi) && y.hi.gte(this.lo);
        }
    };
    /**
     * Return true if the interior of this interval contains any point of the
     * interval 'y' (including its boundary). Works for empty, full, and singleton
     * intervals.
     */
    S1Interval.prototype.interiorIntersects = function (y) {
        if (this.isEmpty() || y.isEmpty() || this.lo.eq(this.hi)) {
            return false;
        }
        if (this.isInverted()) {
            return y.isInverted() || y.lo.lt(this.hi) || y.hi.gt(this.lo);
        }
        else {
            if (y.isInverted()) {
                return y.lo.lt(this.hi) || y.hi.gt(this.lo);
            }
            return (y.lo.lt(this.hi) && y.hi.gt(this.lo)) || this.isFull();
        }
    };
    /**
     * Expand the interval by the minimum amount necessary so that it contains the
     * given point "p" (an angle in the range [-Pi, Pi]).
     */
    S1Interval.prototype.addPoint = function (_p) {
        var p = S2.toDecimal(_p);
        // assert (Math.abs(p) <= S2.M_PI);
        if (p.eq(-S2.M_PI)) {
            p = S2.toDecimal(S2.M_PI);
        }
        if (this.fastContains(p)) {
            return new S1Interval(this.lo, this.hi);
        }
        if (this.isEmpty()) {
            return S1Interval.fromPoint(p);
        }
        else {
            // Compute distance from p to each endpoint.
            var dlo = S1Interval.positiveDistance(p, this.lo);
            var dhi = S1Interval.positiveDistance(this.hi, p);
            if (dlo.lt(dhi)) {
                return new S1Interval(p, this.hi);
            }
            else {
                return new S1Interval(this.lo, p);
            }
            // Adding a point can never turn a non-full interval into a full one.
        }
    };
    /**
     * Return an interval that contains all points within a distance "radius" of
     * a point in this interval. Note that the expansion of an empty interval is
     * always empty. The radius must be non-negative.
     */
    S1Interval.prototype.expanded = function (_radius) {
        var radius = S2.toDecimal(_radius);
        // assert (radius >= 0);
        if (this.isEmpty()) {
            return this;
        }
        // Check whether this interval will be full after expansion, allowing
        // for a 1-bit rounding error when computing each endpoint.
        if (this.getLength().plus(radius.times(2)).gte(2 * S2.M_PI - 1e-15)) {
            return S1Interval.full();
        }
        // NOTE(dbeaumont): Should this remainder be 2 * M_PI or just M_PI ??
        var lo = S2.IEEEremainder(this.lo.minus(radius), 2 * S2.M_PI);
        var hi = S2.IEEEremainder(this.hi.plus(radius), 2 * S2.M_PI);
        if (lo.eq(-S2.M_PI)) {
            lo = S2.toDecimal(S2.M_PI);
        }
        return new S1Interval(lo, hi);
    };
    /**
     * Return the smallest interval that contains this interval and the given
     * interval "y".
     */
    S1Interval.prototype.union = function (y) {
        // The y.is_full() case is handled correctly in all cases by the code
        // below, but can follow three separate code paths depending on whether
        // this interval is inverted, is non-inverted but contains Pi, or neither.
        if (y.isEmpty()) {
            return this;
        }
        if (this.fastContains(y.lo)) {
            if (this.fastContains(y.hi)) {
                // Either this interval contains y, or the union of the two
                // intervals is the Full() interval.
                if (this.containsI(y)) {
                    return this; // is_full() code path
                }
                return S1Interval.full();
            }
            return new S1Interval(this.lo, this.hi, true);
        }
        if (this.fastContains(y.hi)) {
            return new S1Interval(y.lo, this.hi, true);
        }
        // This interval contains neither endpoint of y. This means that either y
        // contains all of this interval, or the two intervals are disjoint.
        if (this.isEmpty() || y.fastContains(this.lo)) {
            return y;
        }
        // Check which pair of endpoints are closer together.
        var dlo = S1Interval.positiveDistance(y.hi, this.lo);
        var dhi = S1Interval.positiveDistance(this.hi, y.lo);
        if (dlo < dhi) {
            return new S1Interval(y.lo, this.hi, true);
        }
        else {
            return new S1Interval(this.lo, y.hi, true);
        }
    };
    /**
     * Return the smallest interval that contains the intersection of this
     * interval with "y". Note that the region of intersection may consist of two
     * disjoint intervals.
     */
    S1Interval.prototype.intersection = function (y) {
        // The y.is_full() case is handled correctly in all cases by the code
        // below, but can follow three separate code paths depending on whether
        // this interval is inverted, is non-inverted but contains Pi, or neither.
        if (y.isEmpty()) {
            return S1Interval.empty();
        }
        if (this.fastContains(y.lo)) {
            if (this.fastContains(y.hi)) {
                // Either this interval contains y, or the region of intersection
                // consists of two disjoint subintervals. In either case, we want
                // to return the shorter of the two original intervals.
                if (y.getLength().lt(this.getLength())) {
                    return y; // is_full() code path
                }
                return this;
            }
            return new S1Interval(y.lo, this.hi, true);
        }
        if (this.fastContains(y.hi)) {
            return new S1Interval(this.lo, y.hi, true);
        }
        // This interval contains neither endpoint of y. This means that either y
        // contains all of this interval, or the two intervals are disjoint.
        if (y.fastContains(this.lo)) {
            return this; // is_empty() okay here
        }
        // assert (!intersects(y));
        return S1Interval.empty();
    };
    /**
     * Return true if the length of the symmetric difference between the two
     * intervals is at most the given tolerance.
     */
    S1Interval.prototype.approxEquals = function (y, maxError) {
        if (maxError === void 0) { maxError = 1e-9; }
        if (this.isEmpty()) {
            return y.getLength().lte(maxError);
        }
        if (y.isEmpty()) {
            return this.getLength().lte(maxError);
        }
        return S2.IEEEremainder(y.lo.minus(this.lo), 2 * S2.M_PI).abs()
            .plus(S2.IEEEremainder(y.hi.minus(this.hi), 2 * S2.M_PI).abs())
            .lte(maxError);
    };
    S1Interval.empty = function () {
        return new S1Interval(S2.M_PI, -S2.M_PI, true);
    };
    S1Interval.full = function () {
        return new S1Interval(-S2.M_PI, S2.M_PI, true);
    };
    S1Interval.fromPoint = function (_p) {
        var p = S2.toDecimal(_p);
        if (p.eq(-S2.M_PI)) {
            p = S2.toDecimal(S2.M_PI);
        }
        return new S1Interval(p, p, true);
    };
    /**
     * Convenience method to construct the minimal interval containing the two
     * given points. This is equivalent to starting with an empty interval and
     * calling AddPoint() twice, but it is more efficient.
     */
    S1Interval.fromPointPair = function (_p1, _p2) {
        // assert (Math.abs(p1) <= S2.M_PI && Math.abs(p2) <= S2.M_PI);
        var p1 = S2.toDecimal(_p1);
        var p2 = S2.toDecimal(_p2);
        if (p1.eq(-S2.M_PI)) {
            p1 = S2.toDecimal(S2.M_PI);
        }
        if (p2.eq(-S2.M_PI)) {
            p2 = S2.toDecimal(S2.M_PI);
        }
        if (S1Interval.positiveDistance(p1, p2).lte(S2.M_PI)) {
            return new S1Interval(p1, p2, true);
        }
        else {
            return new S1Interval(p2, p1, true);
        }
    };
    /**
     * Compute the distance from "a" to "b" in the range [0, 2*Pi). This is
     * equivalent to (drem(b - a - S2.M_PI, 2 * S2.M_PI) + S2.M_PI), except that
     * it is more numerically stable (it does not lose precision for very small
     * positive distances).
     */
    S1Interval.positiveDistance = function (_a, _b) {
        var a = S2.toDecimal(_a);
        var b = S2.toDecimal(_b);
        var d = b.minus(a);
        if (d.gte(0)) {
            return d;
        }
        // We want to ensure that if b == Pi and a == (-Pi + eps),
        // the return result is approximately 2*Pi and not zero.
        return b.plus(S2.M_PI).minus(a.minus(S2.M_PI));
    };
    return S1Interval;
}(Interval));

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
var R1Interval = /** @class */ (function (_super) {
    __extends(R1Interval, _super);
    function R1Interval() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    /** Return true if the interval is empty, i.e. it contains no points. */
    R1Interval.prototype.isEmpty = function () {
        return this.lo.gt(this.hi);
    };
    R1Interval.prototype.getCenter = function () {
        return this.lo.plus(this.hi).dividedBy(2);
    };
    R1Interval.prototype.getLength = function () {
        return this.hi.minus(this.lo);
    };
    R1Interval.prototype.contains = function (_p) {
        var p = S2.toDecimal(_p);
        return p.gte(this.lo) && p.lte(this.hi);
    };
    /** Return true if the interior of the interval contains the point 'p'. */
    R1Interval.prototype.interiorContains = function (_p) {
        var p = S2.toDecimal(_p);
        return p.gt(this.lo) && p.lt(this.hi);
    };
    /**
     * Return true if the interval contains the given interval 'y'. Works for
     * empty, full, and singleton intervals.
     */
    R1Interval.prototype.containsI = function (y) {
        if (y.isEmpty()) {
            return true;
        }
        return y.lo.gte(this.lo) && y.hi.lte(this.hi);
    };
    R1Interval.prototype.interiorContainsI = function (y) {
        if (y.isEmpty()) {
            return true;
        }
        return y.lo.gt(this.lo) && y.hi.lt(this.hi);
    };
    /**
     * Return true if this interval intersects the given interval, i.e. if they
     * have any points in common.
     */
    R1Interval.prototype.intersects = function (y) {
        if (this.lo.lte(y.lo)) {
            return y.lo.lte(this.hi) && y.lo.lte(y.hi);
        }
        else {
            return this.lo.lte(y.hi) && this.lo.lte(this.hi);
        }
    };
    /**
     * Return true if the interior of this interval intersects any point of the
     * given interval (including its boundary).
     */
    R1Interval.prototype.interiorIntersects = function (y) {
        return y.lo.lt(this.hi) && this.lo.lt(y.hi) && this.lo.lt(this.hi) && y.lo.lte(y.hi);
    };
    /** Expand the interval so that it contains the given point "p". */
    R1Interval.prototype.addPoint = function (_p) {
        var p = S2.toDecimal(_p);
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
    };
    /**
     * Return an interval that contains all points with a distance "radius" of a
     * point in this interval. Note that the expansion of an empty interval is
     * always empty.
     */
    R1Interval.prototype.expanded = function (_radius) {
        var radius = S2.toDecimal(_radius);
        // assert (radius >= 0);
        if (this.isEmpty()) {
            return this;
        }
        return new R1Interval(this.lo.minus(radius), this.hi.plus(radius));
    };
    /**
     * Return the smallest interval that contains this interval and the given
     * interval "y".
     */
    R1Interval.prototype.union = function (y) {
        if (this.isEmpty()) {
            return y;
        }
        if (y.isEmpty()) {
            return this;
        }
        return new R1Interval(decimal.Decimal.min(this.lo, y.lo), decimal.Decimal.max(this.hi, y.hi));
    };
    /**
     * Return the intersection of this interval with the given interval. Empty
     * intervals do not need to be special-cased.
     */
    R1Interval.prototype.intersection = function (y) {
        return new R1Interval(decimal.Decimal.max(this.lo, y.lo), decimal.Decimal.min(this.hi, y.hi));
    };
    /**
     * Return true if the length of the symmetric difference between the two
     * intervals is at most the given tolerance.
     */
    R1Interval.prototype.approxEquals = function (y, maxError) {
        if (maxError === void 0) { maxError = 1e-15; }
        if (this.isEmpty()) {
            return y.getLength().lte(maxError);
        }
        if (y.isEmpty()) {
            return this.getLength().lte(maxError);
        }
        return y.lo.minus(this.lo).abs()
            .plus(y.hi.minus(this.hi).abs())
            .lte(maxError);
    };
    R1Interval.empty = function () {
        return new R1Interval(1, 0);
    };
    R1Interval.fromPoint = function (p) {
        return new R1Interval(p, p);
    };
    /**
     * Convenience method to construct the minimal interval containing the two
     * given points. This is equivalent to starting with an empty interval and
     * calling AddPoint() twice, but it is more efficient.
     */
    R1Interval.fromPointPair = function (_p1, _p2) {
        var p1 = S2.toDecimal(_p1);
        var p2 = S2.toDecimal(_p2);
        if (p1.lte(p2)) {
            return new R1Interval(p1, p2);
        }
        else {
            return new R1Interval(p2, p1);
        }
    };
    return R1Interval;
}(Interval));

/*
 * Copyright 2005 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/**
 * This class represents a point on the unit sphere as a pair of
 * latitude-longitude coordinates. Like the rest of the "geometry" package, the
 * intent is to represent spherical geometry as a mathematical abstraction, so
 * functions that are specifically related to the Earth's geometry (e.g.
 * easting/northing conversions) should be put elsewhere.
 *
 */
var S2LatLng = /** @class */ (function () {
    function S2LatLng(latRadians, lngRadians) {
        this.latRadians = S2.toDecimal(latRadians);
        this.lngRadians = S2.toDecimal(lngRadians);
    }
    Object.defineProperty(S2LatLng.prototype, "latDegrees", {
        get: function () {
            return new S1Angle(this.latRadians).degrees();
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(S2LatLng.prototype, "lngDegrees", {
        get: function () {
            return new S1Angle(this.lngRadians).degrees();
        },
        enumerable: true,
        configurable: true
    });
    // Clamps the latitude to the range [-90, 90] degrees, and adds or subtracts
    // a multiple of 360 degrees to the longitude if necessary to reduce it to
    // the range [-180, 180].
    /** Convert an S2LatLng to the equivalent unit-length vector (S2Point). */
    S2LatLng.prototype.toPoint = function () {
        var phi = this.latRadians;
        var theta = this.lngRadians;
        var cosphi = decimal.Decimal.cos(phi);
        return new S2Point(decimal.Decimal.cos(theta).times(cosphi), decimal.Decimal.sin(theta).times(cosphi), decimal.Decimal.sin(phi));
    };
    /**
     * Returns a new S2LatLng based on this instance for which {@link #isValid()}
     * will be {@code true}.
     * <ul>
     * <li>Latitude is clipped to the range {@code [-90, 90]}
     * <li>Longitude is normalized to be in the range {@code [-180, 180]}
     * </ul>
     * <p>If the current point is valid then the returned point will have the same
     * coordinates.
     */
    S2LatLng.prototype.normalized = function () {
        // drem(x, 2 * S2.M_PI) reduces its argument to the range
        // [-S2.M_PI, S2.M_PI] inclusive, which is what we want here.
        return new S2LatLng(decimal.Decimal.max(-S2.M_PI_2, decimal.Decimal.min(S2.M_PI_2, this.latRadians)), S2.IEEEremainder(this.lngRadians, S2.toDecimal(2).times(S2.M_PI)));
        // return new S2LatLng(Math.max(-S2.M_PI_2, Math.min(S2.M_PI_2, this.latRadians)),
        //     S2.IEEEremainder(this.lngRadians, 2 * S2.M_PI));
    };
    S2LatLng.fromDegrees = function (latDegrees, lngDegrees) {
        return new S2LatLng(S1Angle.degrees(latDegrees).radians, S1Angle.degrees(lngDegrees).radians);
    };
    S2LatLng.fromPoint = function (p) {
        return new S2LatLng(S2LatLng.latitude(p).radians, S2LatLng.longitude(p).radians);
    };
    /**
     * Return true if the latitude is between -90 and 90 degrees inclusive and the
     * longitude is between -180 and 180 degrees inclusive.
     */
    S2LatLng.prototype.isValid = function () {
        return this.latRadians.abs().lte(S2.M_PI_2) &&
            this.lngRadians.abs().lte(S2.M_PI);
    };
    /**
     * Scales this point by the given scaling factor.
     * Note that there is no guarantee that the new point will be <em>valid</em>.
     */
    S2LatLng.prototype.mul = function (m) {
        return new S2LatLng(this.latRadians.times(m), this.lngRadians.times(m));
    };
    S2LatLng.latitude = function (p) {
        // We use atan2 rather than asin because the input vector is not necessarily
        // unit length, and atan2 is much more accurate than asin near the poles.
        return new S1Angle(decimal.Decimal.atan2(p.z, p.x.pow(2)
            .plus(p.y.pow(2))
            .sqrt())
        // Math.atan2(p.z, Math.sqrt(p.x * p.x + p.y * p.y))
        );
    };
    S2LatLng.longitude = function (p) {
        // Note that atan2(0, 0) is defined to be zero.
        return new S1Angle(decimal.Decimal.atan2(p.y, p.x));
    };
    S2LatLng.prototype.equals = function (other) {
        return other.latRadians === this.latRadians && other.lngRadians === this.lngRadians;
    };
    S2LatLng.prototype.pointAtDistance = function (_distanceInKm, _bearingRadians) {
        var distanceInM = S2.toDecimal(_distanceInKm).times(1000);
        var distanceToRadius = distanceInM.dividedBy(S2LatLng.EARTH_RADIUS_METERS);
        var bearingRadians = S2.toDecimal(_bearingRadians);
        this.latRadians.sin();
        distanceToRadius.cos();
        var newLat = this.latRadians.sin()
            .times(distanceToRadius.cos())
            .plus(this.latRadians.cos()
            .times(distanceToRadius.sin())
            .times(bearingRadians.cos())).asin();
        var newLng = this.lngRadians
            .plus(decimal.Decimal.atan2(bearingRadians.sin()
            .times(distanceToRadius.sin())
            .times(this.latRadians.cos()), distanceToRadius.cos()
            .minus(this.latRadians.sin().times(newLat.sin()))));
        return new S2LatLng(newLat, newLng);
    };
    /**
     * Generates n LatLngs given a distance in km and the number of points wanted.
     * Generated points will be returned in a Clockwise order starting from North.
     * @param _distanceInKm
     * @param nPoints
     * @returns {S2LatLng[]}
     */
    S2LatLng.prototype.pointsAtDistance = function (_distanceInKm, nPoints) {
        var _this = this;
        if (nPoints === void 0) { nPoints = 4; }
        return Array.apply(null, new Array(nPoints)) // create an array filled of undefined!
            .map(function (p, idx) {
            return S2.toDecimal(360).dividedBy(nPoints).times(idx);
        })
            .map(function (bearingDegree) { return S1Angle.degrees(bearingDegree).radians; })
            .map(function (bearingRadians) { return _this.pointAtDistance(_distanceInKm, bearingRadians); });
    };
    S2LatLng.prototype.getEarthDistance = function (other) {
        return this.getDistance(other).radians.times(S2LatLng.EARTH_RADIUS_METERS);
    };
    S2LatLng.prototype.getDistance = function (other) {
        // This implements the Haversine formula, which is numerically stable for
        // small distances but only gets about 8 digits of precision for very large
        // distances (e.g. antipodal points). Note that 8 digits is still accurate
        // to within about 10cm for a sphere the size of the Earth.
        //
        // This could be fixed with another sin() and cos() below, but at that point
        // you might as well just convert both arguments to S2Points and compute the
        // distance that way (which gives about 15 digits of accuracy for all
        // distances).
        var dLat = other.latRadians.minus(this.latRadians).times(0.5).sin();
        var dLng = other.lngRadians.minus(this.lngRadians).times(0.5).sin();
        var x = dLat.pow(2)
            .plus(dLng.pow(2)
            .times(this.latRadians.cos())
            .times(other.latRadians.cos()));
        // double x = dlat * dlat + dlng * dlng * Math.cos(lat1) * Math.cos(lat2);
        return new S1Angle(S2.toDecimal(2)
            .times(decimal.Decimal.atan2(x.sqrt(), decimal.Decimal.max(0, x.neg().plus(1))
            .sqrt())));
        // Return the distance (measured along the surface of the sphere) to the
        // given S2LatLng. This is mathematically equivalent to:
        //
        // S1Angle::FromRadians(ToPoint().Angle(o.ToPoint())
        //
        // but this implementation is slightly more efficient.
    };
    S2LatLng.prototype.toString = function () {
        return "(" + this.latRadians + ", " + this.lngRadians + ")";
    };
    S2LatLng.prototype.toStringDegrees = function () {
        return "(" + this.latDegrees + ", " + this.lngDegrees + ")";
    };
    S2LatLng.prototype.toGEOJSON = function () {
        return {
            type: 'Feature',
            geometry: {
                type: "Point",
                coordinates: [this.lngDegrees.toNumber(), this.latDegrees.toNumber()]
            },
            properties: {}
        };
    };
    /**
     * Approximate "effective" radius of the Earth in meters.
     */
    S2LatLng.EARTH_RADIUS_METERS = 6367000.0;
    /** The center point the lat/lng coordinate system. */
    S2LatLng.CENTER = new S2LatLng(0.0, 0.0);
    return S2LatLng;
}());

/*
 * Copyright 2006 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/**
 * This class contains various utility functions related to edges. It collects
 * together common code that is needed to implement polygonal geometry such as
 * polylines, loops, and general polygons.
 *
 */
var S2EdgeUtil = /** @class */ (function () {
    function S2EdgeUtil() {
    }
    //   /**
    //    * IEEE floating-point operations have a maximum error of 0.5 ULPS (units in
    //    * the last place). For double-precision numbers, this works out to 2**-53
    //    * (about 1.11e-16) times the magnitude of the result. It is possible to
    //    * analyze the calculation done by getIntersection() and work out the
    //    * worst-case rounding error. I have done a rough version of this, and my
    //    * estimate is that the worst case distance from the intersection point X to
    //    * the great circle through (a0, a1) is about 12 ULPS, or about 1.3e-15. This
    //    * needs to be increased by a factor of (1/0.866) to account for the
    //    * edgeSpliceFraction() in S2PolygonBuilder. Note that the maximum error
    //    * measured by the unittest in 1,000,000 trials is less than 3e-16.
    //    */
    //   public static final S1Angle DEFAULT_INTERSECTION_TOLERANCE = S1Angle.radians(1.5e-15);
    //
    //   /**
    //    * This class allows a vertex chain v0, v1, v2, ... to be efficiently tested
    //    * for intersection with a given fixed edge AB.
    //    */
    //   public static class EdgeCrosser {
    //   // The fields below are all constant.
    //
    //   private final S2Point a;
    //   private final S2Point b;
    //   private final S2Point aCrossB;
    //
    //   // The fields below are updated for each vertex in the chain.
    //
    //   // Previous vertex in the vertex chain.
    //   private S2Point c;
    //   // The orientation of the triangle ACB.
    //   private int acb;
    //
    //   /**
    //    * AB is the given fixed edge, and C is the first vertex of the vertex
    //    * chain. All parameters must point to fixed storage that persists for the
    //    * lifetime of the EdgeCrosser object.
    //    */
    //   public EdgeCrosser(S2Point a, S2Point b, S2Point c) {
    //   this.a = a;
    //   this.b = b;
    //   this.aCrossB = S2Point.crossProd(a, b);
    //   restartAt(c);
    // }
    //
    // /**
    //  * Call this function when your chain 'jumps' to a new place.
    //  */
    // public void restartAt(S2Point c) {
    //   this.c = c;
    //   this.acb = -S2.robustCCW(this.a, this.b, c, this.aCrossB);
    // }
    //
    // /**
    //  * This method is equivalent to calling the S2EdgeUtil.robustCrossing()
    //  * function (defined below) on the edges AB and CD. It returns +1 if there
    //  * is a crossing, -1 if there is no crossing, and 0 if two points from
    //  * different edges are the same. Returns 0 or -1 if either edge is
    //  * degenerate. As a side effect, it saves vertex D to be used as the next
    //  * vertex C.
    //  */
    // public int robustCrossing(S2Point d) {
    //   // For there to be an edge crossing, the triangles ACB, CBD, BDA, DAC must
    //   // all be oriented the same way (CW or CCW). We keep the orientation
    //   // of ACB as part of our state. When each new point D arrives, we
    //   // compute the orientation of BDA and check whether it matches ACB.
    //   // This checks whether the points C and D are on opposite sides of the
    //   // great circle through AB.
    //
    //   // Recall that robustCCW is invariant with respect to rotating its
    //   // arguments, i.e. ABC has the same orientation as BDA.
    //   int bda = S2.robustCCW(this.a, this.b, d, this.aCrossB);
    //   int result;
    //
    //   if (bda == -this.acb && bda != 0) {
    //     // Most common case -- triangles have opposite orientations.
    //     result = -1;
    //   } else if ((bda & this.acb) == 0) {
    //     // At least one value is zero -- two vertices are identical.
    //     result = 0;
    //   } else {
    //     // assert (bda == acb && bda != 0);
    //     result = robustCrossingInternal(d); // Slow path.
    //   }
    //   // Now save the current vertex D as the next vertex C, and also save the
    //   // orientation of the new triangle ACB (which is opposite to the current
    //   // triangle BDA).
    //   this.c = d;
    //   this.acb = -bda;
    //   return result;
    // }
    //
    // /**
    //  * This method is equivalent to the S2EdgeUtil.edgeOrVertexCrossing() method
    //  * defined below. It is similar to robustCrossing, but handles cases where
    //  * two vertices are identical in a way that makes it easy to implement
    //  * point-in-polygon containment tests.
    //  */
    // public boolean edgeOrVertexCrossing(S2Point d) {
    //   // We need to copy c since it is clobbered by robustCrossing().
    //   S2Point c2 = new S2Point(this.c.get(0), this.c.get(1), this.c.get(2));
    //
    //   int crossing = robustCrossing(d);
    //   if (crossing < 0) {
    //     return false;
    //   }
    //   if (crossing > 0) {
    //     return true;
    //   }
    //
    //   return vertexCrossing(this.a, this.b, c2, d);
    // }
    //
    // /**
    //  * This function handles the "slow path" of robustCrossing().
    //  */
    // private int robustCrossingInternal(S2Point d) {
    //   // ACB and BDA have the appropriate orientations, so now we check the
    //   // triangles CBD and DAC.
    //   S2Point cCrossD = S2Point.crossProd(this.c, d);
    //   int cbd = -S2.robustCCW(this.c, d, this.b, cCrossD);
    //   if (cbd != this.acb) {
    //     return -1;
    //   }
    //
    //   int dac = S2.robustCCW(this.c, d, this.a, cCrossD);
    //   return (dac == this.acb) ? 1 : -1;
    // }
    // }
    //
    // /**
    //  * This class computes a bounding rectangle that contains all edges defined by
    //  * a vertex chain v0, v1, v2, ... All vertices must be unit length. Note that
    //  * the bounding rectangle of an edge can be larger than the bounding rectangle
    //  * of its endpoints, e.g. consider an edge that passes through the north pole.
    //  */
    // public static class RectBounder {
    //   // The previous vertex in the chain.
    //   private S2Point a;
    //
    //   // The corresponding latitude-longitude.
    //   private S2LatLng aLatLng;
    //
    //   // The current bounding rectangle.
    //   private S2LatLngRect bound;
    //
    //   public RectBounder() {
    //     this.bound = S2LatLngRect.empty();
    //   }
    //
    //   /**
    //    * This method is called to add each vertex to the chain. 'b' must point to
    //    * fixed storage that persists for the lifetime of the RectBounder.
    //    */
    //   public void addPoint(S2Point b) {
    //   // assert (S2.isUnitLength(b));
    //
    //   S2LatLng bLatLng = new S2LatLng(b);
    //
    //   if (this.bound.isEmpty()) {
    //   this.bound = this.bound.addPoint(bLatLng);
    // } else {
    //   // We can't just call bound.addPoint(bLatLng) here, since we need to
    //   // ensure that all the longitudes between "a" and "b" are included.
    //   this.bound = this.bound.union(S2LatLngRect.fromPointPair(this.aLatLng, bLatLng));
    //
    //   // Check whether the min/max latitude occurs in the edge interior.
    //   // We find the normal to the plane containing AB, and then a vector
    //   // "dir" in this plane that also passes through the equator. We use
    //   // RobustCrossProd to ensure that the edge normal is accurate even
    //   // when the two points are very close together.
    //   S2Point aCrossB = S2.robustCrossProd(this.a, b);
    //   S2Point dir = S2Point.crossProd(aCrossB, new S2Point(0, 0, 1));
    //   double da = dir.dotProd(this.a);
    //   double db = dir.dotProd(b);
    //
    //   if (da * db < 0) {
    //     // Minimum/maximum latitude occurs in the edge interior. This affects
    //     // the latitude bounds but not the longitude bounds.
    //     double absLat = Math.acos(Math.abs(aCrossB.get(2) / aCrossB.norm()));
    //     R1Interval lat = this.bound.lat();
    //     if (da < 0) {
    //       // It's possible that absLat < lat.lo() due to numerical errors.
    //       lat = new R1Interval(lat.lo(), Math.max(absLat, this.bound.lat().hi()));
    //     } else {
    //       lat = new R1Interval(Math.min(-absLat, this.bound.lat().lo()), lat.hi());
    //     }
    //     this.bound = new S2LatLngRect(lat, this.bound.lng());
    //   }
    // }
    // this.a = b;
    // this.aLatLng = bLatLng;
    // }
    //
    // /**
    //  * Return the bounding rectangle of the edge chain that connects the
    //  * vertices defined so far.
    //  */
    // public S2LatLngRect getBound() {
    //   return this.bound;
    // }
    //
    // }
    //
    // /**
    //  * The purpose of this class is to find edges that intersect a given XYZ
    //  * bounding box. It can be used as an efficient rejection test when attempting to
    //  * find edges that intersect a given region. It accepts a vertex chain v0, v1,
    //  * v2, ... and returns a boolean value indicating whether each edge intersects
    //  * the specified bounding box.
    //  *
    //  * We use XYZ intervals instead of something like longitude intervals because
    //  * it is cheap to collect from S2Point lists and any slicing strategy should
    //  * give essentially equivalent results.  See S2Loop for an example of use.
    //  */
    // public static class XYZPruner {
    //   private S2Point lastVertex;
    //
    //   // The region to be tested against.
    //   private boolean boundSet;
    //   private double xmin;
    //   private double ymin;
    //   private double zmin;
    //   private double xmax;
    //   private double ymax;
    //   private double zmax;
    //   private double maxDeformation;
    //
    //   public XYZPruner() {
    //     this.boundSet = false;
    //   }
    //
    //   /**
    //    * Accumulate a bounding rectangle from provided edges.
    //    *
    //    * @param from start of edge
    //    * @param to end of edge.
    //    */
    //   public void addEdgeToBounds(S2Point from, S2Point to) {
    //   if (!this.boundSet) {
    //   this.boundSet = true;
    //   this.xmin = this.xmax = from.x;
    //   this.ymin = this.ymax = from.y;
    //   this.zmin = this.zmax = from.z;
    // }
    // this.xmin = Math.min(this.xmin, Math.min(to.x, from.x));
    // this.ymin = Math.min(this.ymin, Math.min(to.y, from.y));
    // this.zmin = Math.min(this.zmin, Math.min(to.z, from.z));
    // this.xmax = Math.max(this.xmax, Math.max(to.x, from.x));
    // this.ymax = Math.max(this.ymax, Math.max(to.y, from.y));
    // this.zmax = Math.max(this.zmax, Math.max(to.z, from.z));
    //
    // // Because our arcs are really geodesics on the surface of the earth
    // // an edge can have intermediate points outside the xyz bounds implicit
    // // in the end points.  Based on the length of the arc we compute a
    // // generous bound for the maximum amount of deformation.  For small edges
    // // it will be very small but for some large arcs (ie. from (1N,90W) to
    // // (1N,90E) the path can be wildly deformed.  I did a bunch of
    // // experiments with geodesics to get safe bounds for the deformation.
    // double approxArcLen =
    //     Math.abs(from.x - to.x) + Math.abs(from.y - to.y) + Math.abs(from.z - to.z);
    // if (approxArcLen < 0.025) { // less than 2 degrees
    //   this.maxDeformation = Math.max(this.maxDeformation, approxArcLen * 0.0025);
    // } else if (approxArcLen < 1.0) { // less than 90 degrees
    //   this.maxDeformation = Math.max(this.maxDeformation, approxArcLen * 0.11);
    // } else {
    //   this.maxDeformation = approxArcLen * 0.5;
    // }
    // }
    //
    // public void setFirstIntersectPoint(S2Point v0) {
    //   this.xmin = this.xmin - this.maxDeformation;
    //   this.ymin = this.ymin - this.maxDeformation;
    //   this.zmin = this.zmin - this.maxDeformation;
    //   this.xmax = this.xmax + this.maxDeformation;
    //   this.ymax = this.ymax + this.maxDeformation;
    //   this.zmax = this.zmax + this.maxDeformation;
    //   this.lastVertex = v0;
    // }
    //
    // /**
    //  * Returns true if the edge going from the last point to this point passes
    //  * through the pruner bounding box, otherwise returns false.  So the
    //  * method returns false if we are certain there is no intersection, but it
    //  * may return true when there turns out to be no intersection.
    //  */
    // public boolean intersects(S2Point v1) {
    //   boolean result = true;
    //
    //   if ((v1.x < this.xmin && this.lastVertex.x < this.xmin) || (v1.x > this.xmax && this.lastVertex.x > this.xmax)) {
    //     result = false;
    //   } else if ((v1.y < this.ymin && this.lastVertex.y < this.ymin) || (v1.y > this.ymax && this.lastVertex.y > this.ymax)) {
    //     result = false;
    //   } else if ((v1.z < this.zmin && this.lastVertex.z < this.zmin) || (v1.z > this.zmax && this.lastVertex.z > this.zmax)) {
    //     result = false;
    //   }
    //
    //   this.lastVertex = v1;
    //   return result;
    // }
    // }
    //
    // /**
    //  * The purpose of this class is to find edges that intersect a given longitude
    //  * interval. It can be used as an efficient rejection test when attempting to
    //  * find edges that intersect a given region. It accepts a vertex chain v0, v1,
    //  * v2, ... and returns a boolean value indicating whether each edge intersects
    //  * the specified longitude interval.
    //  *
    //  * This class is not currently used as the XYZPruner is preferred for
    //  * S2Loop, but this should be usable in similar circumstances.  Be wary
    //  * of the cost of atan2() in conversions from S2Point to longitude!
    //  */
    // public static class LongitudePruner {
    //   // The interval to be tested against.
    //   private S1Interval interval;
    //
    //   // The longitude of the next v0.
    //   private double lng0;
    //
    //   /**
    //    *'interval' is the longitude interval to be tested against, and 'v0' is
    //    * the first vertex of edge chain.
    //    */
    //   public LongitudePruner(S1Interval interval, S2Point v0) {
    //   this.interval = interval;
    //   this.lng0 = S2LatLng.longitude(v0).radians();
    // }
    //
    // /**
    //  * Returns true if the edge (v0, v1) intersects the given longitude
    //  * interval, and then saves 'v1' to be used as the next 'v0'.
    //  */
    // public boolean intersects(S2Point v1) {
    //   double lng1 = S2LatLng.longitude(v1).radians();
    //   boolean result = this.interval.intersects(S1Interval.fromPointPair(this.lng0, lng1));
    //   this.lng0 = lng1;
    //   return result;
    // }
    // }
    //
    // /**
    //  * A wedge relation's test method accepts two edge chains A=(a0,a1,a2) and
    //  * B=(b0,b1,b2) where a1==b1, and returns either -1, 0, or 1 to indicate the
    //  * relationship between the region to the left of A and the region to the left
    //  * of B. Wedge relations are used to determine the local relationship between
    //  * two polygons that share a common vertex.
    //  *
    //  *  All wedge relations require that a0 != a2 and b0 != b2. Other degenerate
    //  * cases (such as a0 == b2) are handled as expected. The parameter "ab1"
    //  * denotes the common vertex a1 == b1.
    //  */
    // public interface WedgeRelation {
    //   int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2);
    // }
    //
    // public static class WedgeContains implements WedgeRelation {
    //   /**
    //    * Given two edge chains (see WedgeRelation above), this function returns +1
    //    * if the region to the left of A contains the region to the left of B, and
    //    * 0 otherwise.
    //    */
    //   @Override
    //   public int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
    //   // For A to contain B (where each loop interior is defined to be its left
    //   // side), the CCW edge order around ab1 must be a2 b2 b0 a0. We split
    //   // this test into two parts that test three vertices each.
    //   return S2.orderedCCW(a2, b2, b0, ab1) && S2.orderedCCW(b0, a0, a2, ab1) ? 1 : 0;
    // }
    // }
    //
    // public static class WedgeIntersects implements WedgeRelation {
    //   /**
    //    * Given two edge chains (see WedgeRelation above), this function returns -1
    //    * if the region to the left of A intersects the region to the left of B,
    //    * and 0 otherwise. Note that regions are defined such that points along a
    //    * boundary are contained by one side or the other, not both. So for
    //    * example, if A,B,C are distinct points ordered CCW around a vertex O, then
    //    * the wedges BOA, AOC, and COB do not intersect.
    //    */
    //   @Override
    //   public int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
    //   // For A not to intersect B (where each loop interior is defined to be
    //   // its left side), the CCW edge order around ab1 must be a0 b2 b0 a2.
    //   // Note that it's important to write these conditions as negatives
    //   // (!OrderedCCW(a,b,c,o) rather than Ordered(c,b,a,o)) to get correct
    //   // results when two vertices are the same.
    //   return (S2.orderedCCW(a0, b2, b0, ab1) && S2.orderedCCW(b0, a2, a0, ab1) ? 0 : -1);
    // }
    // }
    //
    // public static class WedgeContainsOrIntersects implements WedgeRelation {
    //   /**
    //    * Given two edge chains (see WedgeRelation above), this function returns +1
    //    * if A contains B, 0 if A and B are disjoint, and -1 if A intersects but
    //    * does not contain B.
    //    */
    //   @Override
    //   public int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
    //   // This is similar to WedgeContainsOrCrosses, except that we want to
    //   // distinguish cases (1) [A contains B], (3) [A and B are disjoint],
    //   // and (2,4,5,6) [A intersects but does not contain B].
    //
    //   if (S2.orderedCCW(a0, a2, b2, ab1)) {
    //   // We are in case 1, 5, or 6, or case 2 if a2 == b2.
    //   return S2.orderedCCW(b2, b0, a0, ab1) ? 1 : -1; // Case 1 vs. 2,5,6.
    // }
    // // We are in cases 2, 3, or 4.
    // if (!S2.orderedCCW(a2, b0, b2, ab1)) {
    //   return 0; // Case 3.
    // }
    //
    // // We are in case 2 or 4, or case 3 if a2 == b0.
    // return (a2.equals(b0)) ? 0 : -1; // Case 3 vs. 2,4.
    // }
    // }
    //
    // public static class WedgeContainsOrCrosses implements WedgeRelation {
    //   /**
    //    * Given two edge chains (see WedgeRelation above), this function returns +1
    //    * if A contains B, 0 if B contains A or the two wedges do not intersect,
    //    * and -1 if the edge chains A and B cross each other (i.e. if A intersects
    //    * both the interior and exterior of the region to the left of B). In
    //    * degenerate cases where more than one of these conditions is satisfied,
    //    * the maximum possible result is returned. For example, if A == B then the
    //    * result is +1.
    //    */
    //   @Override
    //   public int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
    //   // There are 6 possible edge orderings at a shared vertex (all
    //   // of these orderings are circular, i.e. abcd == bcda):
    //   //
    //   // (1) a2 b2 b0 a0: A contains B
    //   // (2) a2 a0 b0 b2: B contains A
    //   // (3) a2 a0 b2 b0: A and B are disjoint
    //   // (4) a2 b0 a0 b2: A and B intersect in one wedge
    //   // (5) a2 b2 a0 b0: A and B intersect in one wedge
    //   // (6) a2 b0 b2 a0: A and B intersect in two wedges
    //   //
    //   // In cases (4-6), the boundaries of A and B cross (i.e. the boundary
    //   // of A intersects the interior and exterior of B and vice versa).
    //   // Thus we want to distinguish cases (1), (2-3), and (4-6).
    //   //
    //   // Note that the vertices may satisfy more than one of the edge
    //   // orderings above if two or more vertices are the same. The tests
    //   // below are written so that we take the most favorable
    //   // interpretation, i.e. preferring (1) over (2-3) over (4-6). In
    //   // particular note that if orderedCCW(a,b,c,o) returns true, it may be
    //   // possible that orderedCCW(c,b,a,o) is also true (if a == b or b == c).
    //
    //   if (S2.orderedCCW(a0, a2, b2, ab1)) {
    //   // The cases with this vertex ordering are 1, 5, and 6,
    //   // although case 2 is also possible if a2 == b2.
    //   if (S2.orderedCCW(b2, b0, a0, ab1)) {
    //   return 1; // Case 1 (A contains B)
    // }
    //
    // // We are in case 5 or 6, or case 2 if a2 == b2.
    // return (a2.equals(b2)) ? 0 : -1; // Case 2 vs. 5,6.
    // }
    // // We are in case 2, 3, or 4.
    // return S2.orderedCCW(a0, b0, a2, ab1) ? 0 : -1; // Case 2,3 vs. 4.
    // }
    // }
    //
    // /**
    //  * Return true if edge AB crosses CD at a point that is interior to both
    //  * edges. Properties:
    //  *
    //  *  (1) simpleCrossing(b,a,c,d) == simpleCrossing(a,b,c,d) (2)
    //  * simpleCrossing(c,d,a,b) == simpleCrossing(a,b,c,d)
    //  */
    // public static boolean simpleCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    //   // We compute simpleCCW() for triangles ACB, CBD, BDA, and DAC. All
    //   // of these triangles need to have the same orientation (CW or CCW)
    //   // for an intersection to exist. Note that this is slightly more
    //   // restrictive than the corresponding definition for planar edges,
    //   // since we need to exclude pairs of line segments that would
    //   // otherwise "intersect" by crossing two antipodal points.
    //
    //   S2Point ab = S2Point.crossProd(a, b);
    //   double acb = -(ab.dotProd(c));
    //   double bda = ab.dotProd(d);
    //   if (acb * bda <= 0) {
    //     return false;
    //   }
    //
    //   S2Point cd = S2Point.crossProd(c, d);
    //   double cbd = -(cd.dotProd(b));
    //   double dac = cd.dotProd(a);
    //   return (acb * cbd > 0) && (acb * dac > 0);
    // }
    //
    // /**
    //  * Like SimpleCrossing, except that points that lie exactly on a line are
    //  * arbitrarily classified as being on one side or the other (according to the
    //  * rules of S2.robustCCW). It returns +1 if there is a crossing, -1 if there
    //  * is no crossing, and 0 if any two vertices from different edges are the
    //  * same. Returns 0 or -1 if either edge is degenerate. Properties of
    //  * robustCrossing:
    //  *
    //  *  (1) robustCrossing(b,a,c,d) == robustCrossing(a,b,c,d) (2)
    //  * robustCrossing(c,d,a,b) == robustCrossing(a,b,c,d) (3)
    //  * robustCrossing(a,b,c,d) == 0 if a==c, a==d, b==c, b==d (3)
    //  * robustCrossing(a,b,c,d) <= 0 if a==b or c==d
    //  *
    //  *  Note that if you want to check an edge against a *chain* of other edges,
    //  * it is much more efficient to use an EdgeCrosser (above).
    //  */
    // public static int robustCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    //   // For there to be a crossing, the triangles ACB, CBD, BDA, DAC must
    //   // all have the same orientation (clockwise or counterclockwise).
    //   //
    //   // First we compute the orientation of ACB and BDA. We permute the
    //   // arguments to robustCCW so that we can reuse the cross-product of A and B.
    //   // Recall that when the arguments to robustCCW are permuted, the sign of the
    //   // result changes according to the sign of the permutation. Thus ACB and
    //   // ABC are oppositely oriented, while BDA and ABD are the same.
    //   S2Point aCrossB = S2Point.crossProd(a, b);
    //   int acb = -S2.robustCCW(a, b, c, aCrossB);
    //   int bda = S2.robustCCW(a, b, d, aCrossB);
    //
    //   // If any two vertices are the same, the result is degenerate.
    //   if ((bda & acb) == 0) {
    //     return 0;
    //   }
    //
    //   // If ABC and BDA have opposite orientations (the most common case),
    //   // there is no crossing.
    //   if (bda != acb) {
    //     return -1;
    //   }
    //
    //   // Otherwise we compute the orientations of CBD and DAC, and check whether
    //   // their orientations are compatible with the other two triangles.
    //   S2Point cCrossD = S2Point.crossProd(c, d);
    //   int cbd = -S2.robustCCW(c, d, b, cCrossD);
    //   if (cbd != acb) {
    //     return -1;
    //   }
    //
    //   int dac = S2.robustCCW(c, d, a, cCrossD);
    //   return (dac == acb) ? 1 : -1;
    // }
    //
    // /**
    //  * Given two edges AB and CD where at least two vertices are identical (i.e.
    //  * robustCrossing(a,b,c,d) == 0), this function defines whether the two edges
    //  * "cross" in a such a way that point-in-polygon containment tests can be
    //  * implemented by counting the number of edge crossings. The basic rule is
    //  * that a "crossing" occurs if AB is encountered after CD during a CCW sweep
    //  * around the shared vertex starting from a fixed reference point.
    //  *
    //  *  Note that according to this rule, if AB crosses CD then in general CD does
    //  * not cross AB. However, this leads to the correct result when counting
    //  * polygon edge crossings. For example, suppose that A,B,C are three
    //  * consecutive vertices of a CCW polygon. If we now consider the edge
    //  * crossings of a segment BP as P sweeps around B, the crossing number changes
    //  * parity exactly when BP crosses BA or BC.
    //  *
    //  *  Useful properties of VertexCrossing (VC):
    //  *
    //  *  (1) VC(a,a,c,d) == VC(a,b,c,c) == false (2) VC(a,b,a,b) == VC(a,b,b,a) ==
    //  * true (3) VC(a,b,c,d) == VC(a,b,d,c) == VC(b,a,c,d) == VC(b,a,d,c) (3) If
    //  * exactly one of a,b equals one of c,d, then exactly one of VC(a,b,c,d) and
    //  * VC(c,d,a,b) is true
    //  *
    //  * It is an error to call this method with 4 distinct vertices.
    //  */
    // public static boolean vertexCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    //   // If A == B or C == D there is no intersection. We need to check this
    //   // case first in case 3 or more input points are identical.
    //   if (a.equals(b) || c.equals(d)) {
    //     return false;
    //   }
    //
    //   // If any other pair of vertices is equal, there is a crossing if and only
    //   // if orderedCCW() indicates that the edge AB is further CCW around the
    //   // shared vertex than the edge CD.
    //   if (a.equals(d)) {
    //     return S2.orderedCCW(S2.ortho(a), c, b, a);
    //   }
    //   if (b.equals(c)) {
    //     return S2.orderedCCW(S2.ortho(b), d, a, b);
    //   }
    //   if (a.equals(c)) {
    //     return S2.orderedCCW(S2.ortho(a), d, b, a);
    //   }
    //   if (b.equals(d)) {
    //     return S2.orderedCCW(S2.ortho(b), c, a, b);
    //   }
    //
    //   // assert (false);
    //   return false;
    // }
    //
    // /**
    //  * A convenience function that calls robustCrossing() to handle cases where
    //  * all four vertices are distinct, and VertexCrossing() to handle cases where
    //  * two or more vertices are the same. This defines a crossing function such
    //  * that point-in-polygon containment tests can be implemented by simply
    //  * counting edge crossings.
    //  */
    // public static boolean edgeOrVertexCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    //   int crossing = robustCrossing(a, b, c, d);
    //   if (crossing < 0) {
    //     return false;
    //   }
    //   if (crossing > 0) {
    //     return true;
    //   }
    //   return vertexCrossing(a, b, c, d);
    // }
    //
    // static class CloserResult {
    //   private double dmin2;
    //   private S2Point vmin;
    //
    //   public double getDmin2() {
    //   return this.dmin2;
    // }
    //
    //   public S2Point getVmin() {
    //   return this.vmin;
    // }
    //
    //   public CloserResult(double dmin2, S2Point vmin) {
    //   this.dmin2 = dmin2;
    //   this.vmin = vmin;
    // }
    //
    // public void replaceIfCloser(S2Point x, S2Point y) {
    //   // If the squared distance from x to y is less than dmin2, then replace
    //   // vmin by y and update dmin2 accordingly.
    //   double d2 = S2Point.minus(x, y).norm2();
    //   if (d2 < this.dmin2 || (d2 == this.dmin2 && y.lessThan(this.vmin))) {
    //     this.dmin2 = d2;
    //     this.vmin = y;
    //   }
    // }
    // }
    //
    // /*
    //  * Given two edges AB and CD such that robustCrossing() is true, return their
    //  * intersection point. Useful properties of getIntersection (GI):
    //  *
    //  * (1) GI(b,a,c,d) == GI(a,b,d,c) == GI(a,b,c,d) (2) GI(c,d,a,b) ==
    //  * GI(a,b,c,d)
    //  *
    //  * The returned intersection point X is guaranteed to be close to the edges AB
    //  * and CD, but if the edges intersect at a very small angle then X may not be
    //  * close to the true mathematical intersection point P. See the description of
    //  * "DEFAULT_INTERSECTION_TOLERANCE" below for details.
    //  */
    // public static S2Point getIntersection(S2Point a0, S2Point a1, S2Point b0, S2Point b1) {
    //   Preconditions.checkArgument(robustCrossing(a0, a1, b0, b1) > 0,
    //       "Input edges a0a1 and b0b1 muct have a true robustCrossing.");
    //
    //   // We use robustCrossProd() to get accurate results even when two endpoints
    //   // are close together, or when the two line segments are nearly parallel.
    //   S2Point aNorm = S2Point.normalize(S2.robustCrossProd(a0, a1));
    //   S2Point bNorm = S2Point.normalize(S2.robustCrossProd(b0, b1));
    //   S2Point x = S2Point.normalize(S2.robustCrossProd(aNorm, bNorm));
    //
    //   // Make sure the intersection point is on the correct side of the sphere.
    //   // Since all vertices are unit length, and edges are less than 180 degrees,
    //   // (a0 + a1) and (b0 + b1) both have positive dot product with the
    //   // intersection point. We use the sum of all vertices to make sure that the
    //   // result is unchanged when the edges are reversed or exchanged.
    //   if (x.dotProd(S2Point.add(S2Point.add(a0, a1), S2Point.add(b0, b1))) < 0) {
    //     x = S2Point.neg(x);
    //   }
    //
    //   // The calculation above is sufficient to ensure that "x" is within
    //   // DEFAULT_INTERSECTION_TOLERANCE of the great circles through (a0,a1) and
    //   // (b0,b1).
    //   // However, if these two great circles are very close to parallel, it is
    //   // possible that "x" does not lie between the endpoints of the given line
    //   // segments. In other words, "x" might be on the great circle through
    //   // (a0,a1) but outside the range covered by (a0,a1). In this case we do
    //   // additional clipping to ensure that it does.
    //
    //   if (S2.orderedCCW(a0, x, a1, aNorm) && S2.orderedCCW(b0, x, b1, bNorm)) {
    //     return x;
    //   }
    //
    //   // Find the acceptable endpoint closest to x and return it. An endpoint is
    //   // acceptable if it lies between the endpoints of the other line segment.
    //   CloserResult r = new CloserResult(10, x);
    //   if (S2.orderedCCW(b0, a0, b1, bNorm)) {
    //     r.replaceIfCloser(x, a0);
    //   }
    //   if (S2.orderedCCW(b0, a1, b1, bNorm)) {
    //     r.replaceIfCloser(x, a1);
    //   }
    //   if (S2.orderedCCW(a0, b0, a1, aNorm)) {
    //     r.replaceIfCloser(x, b0);
    //   }
    //   if (S2.orderedCCW(a0, b1, a1, aNorm)) {
    //     r.replaceIfCloser(x, b1);
    //   }
    //   return r.getVmin();
    // }
    //
    // /**
    //  * Given a point X and an edge AB, return the distance ratio AX / (AX + BX).
    //  * If X happens to be on the line segment AB, this is the fraction "t" such
    //  * that X == Interpolate(A, B, t). Requires that A and B are distinct.
    //  */
    // public static double getDistanceFraction(S2Point x, S2Point a0, S2Point a1) {
    //   Preconditions.checkArgument(!a0.equals(a1));
    //   double d0 = x.angle(a0);
    //   double d1 = x.angle(a1);
    //   return d0 / (d0 + d1);
    // }
    //
    // /**
    //  * Return the minimum distance from X to any point on the edge AB. The result
    //  * is very accurate for small distances but may have some numerical error if
    //  * the distance is large (approximately Pi/2 or greater). The case A == B is
    //  * handled correctly. Note: x, a and b must be of unit length. Throws
    //  * IllegalArgumentException if this is not the case.
    //  */
    // public static getDistance(x:S2Point , a:S2Point , b:S2Point ):S1Angle  {
    //   return this.getDistance(x, a, b, S2.robustCrossProd(a, b));
    // }
    /**
     * A slightly more efficient version of getDistance() where the cross product
     * of the two endpoints has been precomputed. The cross product does not need
     * to be normalized, but should be computed using S2.robustCrossProd() for the
     * most accurate results.
     */
    S2EdgeUtil.getDistance = function (x, a, b, aCrossB) {
        // Preconditions.checkArgument(S2.isUnitLength(x));
        // Preconditions.checkArgument(S2.isUnitLength(a));
        // Preconditions.checkArgument(S2.isUnitLength(b));
        if (aCrossB === void 0) { aCrossB = S2.robustCrossProd(a, b); }
        // There are three cases. If X is located in the spherical wedge defined by
        // A, B, and the axis A x B, then the closest point is on the segment AB.
        // Otherwise the closest point is either A or B; the dividing line between
        // these two cases is the great circle passing through (A x B) and the
        // midpoint of AB.
        if (S2.simpleCCW(aCrossB, a, x) && S2.simpleCCW(x, b, aCrossB)) {
            // The closest point to X lies on the segment AB. We compute the distance
            // to the corresponding great circle. The result is accurate for small
            // distances but not necessarily for large distances (approaching Pi/2).
            var sinDist = x.dotProd(aCrossB).abs().dividedBy(aCrossB.norm());
            return new S1Angle(decimal.Decimal.asin(decimal.Decimal.min(1.0, sinDist)));
        }
        // Otherwise, the closest point is either A or B. The cheapest method is
        // just to compute the minimum of the two linear (as opposed to spherical)
        // distances and convert the result to an angle. Again, this method is
        // accurate for small but not large distances (approaching Pi).
        var linearDist2 = decimal.Decimal.min(S2Point.minus(x, a).norm2(), S2Point.minus(x, b).norm2());
        return new S1Angle(decimal.Decimal.asin(decimal.Decimal.min(1.0, linearDist2.sqrt().times(0.5))).times(2));
    };
    return S2EdgeUtil;
}());

var S2LatLngRect = /** @class */ (function () {
    function S2LatLngRect(lat, lng) {
        this.lat = lat;
        this.lng = lng;
    }
    S2LatLngRect.fromLatLng = function (lo, hi) {
        return new S2LatLngRect(new R1Interval(lo.latRadians, hi.latRadians), new S1Interval(lo.lngRadians, hi.lngRadians));
    };
    /** The canonical empty rectangle */
    S2LatLngRect.empty = function () {
        return new S2LatLngRect(R1Interval.empty(), S1Interval.empty());
    };
    /** The canonical full rectangle. */
    S2LatLngRect.full = function () {
        return new S2LatLngRect(S2LatLngRect.fullLat(), S1Interval.full());
    };
    /** The full allowable range of latitudes. */
    S2LatLngRect.fullLat = function () {
        return new R1Interval(-S2.M_PI_2, S2.M_PI_2);
    };
    /**
     * Construct a rectangle from a center point (in lat-lng space) and size in
     * each dimension. If size.lng is greater than 360 degrees it is clamped,
     * and latitudes greater than +/- 90 degrees are also clamped. So for example,
     * FromCenterSize((80,170),(20,20)) -> (lo=(60,150),hi=(90,-170)).
     */
    S2LatLngRect.fromCenterSize = function (center, size) {
        return S2LatLngRect.fromPoint(center).expanded(size.mul(0.5));
    };
    /** Convenience method to construct a rectangle containing a single point. */
    S2LatLngRect.fromPoint = function (p) {
        // assert (p.isValid());
        return S2LatLngRect.fromLatLng(p, p);
    };
    /**
     * Convenience method to construct the minimal bounding rectangle containing
     * the two given points. This is equivalent to starting with an empty
     * rectangle and calling AddPoint() twice. Note that it is different than the
     * S2LatLngRect(lo, hi) constructor, where the first point is always used as
     * the lower-left corner of the resulting rectangle.
     */
    S2LatLngRect.fromPointPair = function (p1, p2) {
        // assert (p1.isValid() && p2.isValid());
        return new S2LatLngRect(R1Interval.fromPointPair(p1.latRadians, p2
            .latRadians), S1Interval.fromPointPair(p1.lngRadians, p2.lngRadians));
    };
    /**
     * Return a latitude-longitude rectangle that contains the edge from "a" to
     * "b". Both points must be unit-length. Note that the bounding rectangle of
     * an edge can be larger than the bounding rectangle of its endpoints.
     */
    S2LatLngRect.fromEdge = function (a, b) {
        // assert (S2.isUnitLength(a) && S2.isUnitLength(b));
        var r = S2LatLngRect.fromPointPair(S2LatLng.fromPoint(a), S2LatLng.fromPoint(b));
        // Check whether the min/max latitude occurs in the edge interior.
        // We find the normal to the plane containing AB, and then a vector "dir" in
        // this plane that also passes through the equator. We use RobustCrossProd
        // to ensure that the edge normal is accurate even when the two points are
        // very close together.
        var ab = S2.robustCrossProd(a, b);
        var dir = S2Point.crossProd(ab, new S2Point(0, 0, 1));
        var da = dir.dotProd(a);
        var db = dir.dotProd(b);
        if (da.times(db).gte(0)) {
            // Minimum and maximum latitude are attained at the vertices.
            return r;
        }
        // Minimum/maximum latitude occurs in the edge interior. This affects the
        // latitude bounds but not the longitude bounds.
        var absLat = decimal.Decimal.acos(ab.z.dividedBy(ab.norm()).abs());
        if (da.lt(0)) {
            return new S2LatLngRect(new R1Interval(r.lat.lo, absLat), r.lng);
        }
        else {
            return new S2LatLngRect(new R1Interval(-absLat, r.lat.hi), r.lng);
        }
    };
    /**
     * Return true if the rectangle is valid, which essentially just means that
     * the latitude bounds do not exceed Pi/2 in absolute value and the longitude
     * bounds do not exceed Pi in absolute value.
     *
     */
    S2LatLngRect.prototype.isValid = function () {
        // The lat/lng ranges must either be both empty or both non-empty.
        return (this.lat.lo.abs().lte(S2.M_PI_2) && this.lat.hi.abs().lte(S2.M_PI_2)
            && this.lng.isValid() && this.lat.isEmpty() == this.lng.isEmpty());
    };
    S2LatLngRect.prototype.lo = function () {
        return new S2LatLng(this.lat.lo, this.lng.lo);
    };
    S2LatLngRect.prototype.hi = function () {
        return new S2LatLng(this.lat.hi, this.lng.hi);
    };
    /**
     * Return true if the rectangle is empty, i.e. it contains no points at all.
     */
    S2LatLngRect.prototype.isEmpty = function () {
        return this.lat.isEmpty();
    };
    // Return true if the rectangle is full, i.e. it contains all points.
    S2LatLngRect.prototype.isFull = function () {
        // console.log(this.lat.toString());
        // console.log(S2LatLngRect.fullLat().toString());
        return this.lat.equals(S2LatLngRect.fullLat()) && this.lng.isFull();
    };
    /**
     * Return true if lng_.lo() > lng_.hi(), i.e. the rectangle crosses the 180
     * degree latitude line.
     */
    S2LatLngRect.prototype.isInverted = function () {
        return this.lng.isInverted();
    };
    /** Return the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order. */
    S2LatLngRect.prototype.getVertex = function (k) {
        // Return the points in CCW order (SW, SE, NE, NW).
        switch (k) {
            case 0:
                return this.lo();
            case 1:
                return new S2LatLng(this.lat.lo, this.lng.hi);
            case 2:
                return this.hi();
            case 3:
                return new S2LatLng(this.lat.hi, this.lng.lo);
            default:
                throw new Error("Invalid vertex index.");
        }
    };
    /**
     * Return the center of the rectangle in latitude-longitude space (in general
     * this is not the center of the region on the sphere).
     */
    S2LatLngRect.prototype.getCenter = function () {
        return new S2LatLng(this.lat.getCenter(), this.lng.getCenter());
    };
    /**
     * Return the minimum distance (measured along the surface of the sphere)
     * from a given point to the rectangle (both its boundary and its interior).
     * The latLng must be valid.
     */
    S2LatLngRect.prototype.getDistanceLL = function (p) {
        // The algorithm here is the same as in getDistance(S2LagLngRect), only
        // with simplified calculations.
        var a = this;
        if (a.isEmpty()) {
            throw new Error();
        }
        if (!p.isValid()) {
            throw new Error('point is not valid');
        }
        if (a.lng.contains(p.lngRadians)) {
            return new S1Angle(decimal.Decimal.max(0.0, decimal.Decimal.max(p.latRadians.minus(a.lat.hi), a.lat.lo.minus(p.latRadians))));
        }
        var interval = new S1Interval(a.lng.hi, a.lng.complement().getCenter());
        var aLng = a.lng.lo;
        if (interval.contains(p.lngRadians)) {
            aLng = a.lng.hi;
        }
        var lo = new S2LatLng(a.lat.lo, aLng).toPoint();
        var hi = new S2LatLng(a.lat.hi, aLng).toPoint();
        var loCrossHi = new S2LatLng(0, aLng.minus(S2.M_PI_2)).normalized().toPoint();
        return S2EdgeUtil.getDistance(p.toPoint(), lo, hi, loCrossHi);
    };
    /**
     * Return the minimum distance (measured along the surface of the sphere) to
     * the given S2LatLngRect. Both S2LatLngRects must be non-empty.
     */
    S2LatLngRect.prototype.getDistanceLLR = function (other) {
        var a = this;
        var b = other;
        if (a.isEmpty()) {
            throw new Error();
        }
        if (b.isEmpty()) {
            throw new Error();
        }
        // First, handle the trivial cases where the longitude intervals overlap.
        if (a.lng.intersects(b.lng)) {
            if (a.lat.intersects(b.lat)) {
                return new S1Angle(0); // Intersection between a and b.
            }
            // We found an overlap in the longitude interval, but not in the latitude
            // interval. This means the shortest path travels along some line of
            // longitude connecting the high-latitude of the lower rect with the
            // low-latitude of the higher rect.
            var lo = void 0, hi = void 0;
            if (a.lat.lo.gt(b.lat.hi)) {
                lo = b.lat.hi;
                hi = a.lat.lo;
            }
            else {
                lo = a.lat.hi;
                hi = b.lat.lo;
            }
            return new S1Angle(hi.radians().minus(lo.radians()));
        }
        // The longitude intervals don't overlap. In this case, the closest points
        // occur somewhere on the pair of longitudinal edges which are nearest in
        // longitude-space.
        var aLng, bLng;
        var loHi = S1Interval.fromPointPair(a.lng.lo, b.lng.hi);
        var hiLo = S1Interval.fromPointPair(a.lng.hi, b.lng.lo);
        if (loHi.getLength().lt(hiLo.getLength())) {
            aLng = a.lng.lo;
            bLng = b.lng.hi;
        }
        else {
            aLng = a.lng.hi;
            bLng = b.lng.lo;
        }
        // The shortest distance between the two longitudinal segments will include
        // at least one segment endpoint. We could probably narrow this down further
        // to a single point-edge distance by comparing the relative latitudes of the
        // endpoints, but for the sake of clarity, we'll do all four point-edge
        // distance tests.
        var aLo = new S2LatLng(a.lat.lo, aLng).toPoint();
        var aHi = new S2LatLng(a.lat.hi, aLng).toPoint();
        var aLoCrossHi = new S2LatLng(0, aLng.radians().minus(S2.M_PI_2)).normalized().toPoint();
        var bLo = new S2LatLng(b.lat.lo, bLng).toPoint();
        var bHi = new S2LatLng(b.lat.hi, bLng).toPoint();
        var bLoCrossHi = new S2LatLng(0, bLng.radians().minus(S2.M_PI_2)).normalized().toPoint();
        return S1Angle.min(S2EdgeUtil.getDistance(aLo, bLo, bHi, bLoCrossHi), S1Angle.min(S2EdgeUtil.getDistance(aHi, bLo, bHi, bLoCrossHi), S1Angle.min(S2EdgeUtil.getDistance(bLo, aLo, aHi, aLoCrossHi), S2EdgeUtil.getDistance(bHi, aLo, aHi, aLoCrossHi))));
    };
    /**
     * Return the width and height of this rectangle in latitude-longitude space.
     * Empty rectangles have a negative width and height.
     */
    S2LatLngRect.prototype.getSize = function () {
        return new S2LatLng(this.lat.getLength(), this.lng.getLength());
    };
    /**
     * More efficient version of Contains() that accepts a S2LatLng rather than an
     * S2Point.
     */
    S2LatLngRect.prototype.containsLL = function (ll) {
        // assert (ll.isValid());
        return (this.lat.contains(ll.latRadians) && this.lng.contains(ll.lngRadians));
    };
    /**
     * Return true if and only if the given point is contained in the interior of
     * the region (i.e. the region excluding its boundary). The point 'p' does not
     * need to be normalized.
     */
    S2LatLngRect.prototype.interiorContainsP = function (p) {
        return this.interiorContainsLL(S2LatLng.fromPoint(p));
    };
    /**
     * More efficient version of InteriorContains() that accepts a S2LatLng rather
     * than an S2Point.
     */
    S2LatLngRect.prototype.interiorContainsLL = function (ll) {
        // assert (ll.isValid());
        return (this.lat.interiorContains(ll.latRadians) && this.lng
            .interiorContains(ll.lngRadians));
    };
    /**
     * Return true if and only if the rectangle contains the given other
     * rectangle.
     */
    S2LatLngRect.prototype.containsLLR = function (other) {
        return this.lat.containsI(other.lat) && this.lng.containsI(other.lng);
    };
    /**
     * Return true if and only if the interior of this rectangle contains all
     * points of the given other rectangle (including its boundary).
     */
    S2LatLngRect.prototype.interiorContainsLLR = function (other) {
        return (this.lat.interiorContainsI(other.lat) && this.lng
            .interiorContainsI(other.lng));
    };
    /** Return true if this rectangle and the given other rectangle have any
     points in common. */
    S2LatLngRect.prototype.intersectsLLR = function (other) {
        return this.lat.intersects(other.lat) && this.lng.intersects(other.lng);
    };
    /**
     * Returns true if this rectangle intersects the given cell. (This is an exact
     * test and may be fairly expensive, see also MayIntersect below.)
     */
    S2LatLngRect.prototype.intersects = function (cell) {
        // First we eliminate the cases where one region completely contains the
        // other. Once these are disposed of, then the regions will intersect
        // if and only if their boundaries intersect.
        if (this.isEmpty()) {
            return false;
        }
        if (this.containsP(cell.getCenter())) {
            return true;
        }
        if (cell.contains(this.getCenter().toPoint())) {
            return true;
        }
        // Quick rejection test (not required for correctness).
        if (!this.intersectsLLR(cell.getRectBound())) {
            return false;
        }
        // Now check whether the boundaries intersect. Unfortunately, a
        // latitude-longitude rectangle does not have straight edges -- two edges
        // are curved, and at least one of them is concave.
        // Precompute the cell vertices as points and latitude-longitudes.
        var cellV = new Array(4);
        var cellLl = new Array(4);
        for (var i = 0; i < 4; ++i) {
            cellV[i] = cell.getVertex(i); // Must be normalized.
            cellLl[i] = S2LatLng.fromPoint(cellV[i]);
            if (this.containsLL(cellLl[i])) {
                return true; // Quick acceptance test.
            }
        }
        for (var i = 0; i < 4; ++i) {
            var edgeLng = S1Interval.fromPointPair(cellLl[i].lngRadians, cellLl[(i + 1) & 3].lngRadians);
            if (!this.lng.intersects(edgeLng)) {
                continue;
            }
            var a = cellV[i];
            var b = cellV[(i + 1) & 3];
            if (edgeLng.contains(this.lng.lo)) {
                if (S2LatLngRect.intersectsLngEdge(a, b, this.lat, this.lng.lo)) {
                    return true;
                }
            }
            if (edgeLng.contains(this.lng.hi)) {
                if (S2LatLngRect.intersectsLngEdge(a, b, this.lat, this.lng.hi)) {
                    return true;
                }
            }
            if (S2LatLngRect.intersectsLatEdge(a, b, this.lat.lo, this.lng)) {
                return true;
            }
            if (S2LatLngRect.intersectsLatEdge(a, b, this.lat.hi, this.lng)) {
                return true;
            }
        }
        return false;
    };
    /**
     * Return true if and only if the interior of this rectangle intersects any
     * point (including the boundary) of the given other rectangle.
     */
    S2LatLngRect.prototype.interiorIntersects = function (other) {
        return (this.lat.interiorIntersects(other.lat) && this.lng
            .interiorIntersects(other.lng));
    };
    S2LatLngRect.prototype.addPoint = function (p) {
        return this.addPointLL(S2LatLng.fromPoint(p));
    };
    // Increase the size of the bounding rectangle to include the given point.
    // The rectangle is expanded by the minimum amount possible.
    S2LatLngRect.prototype.addPointLL = function (ll) {
        var newLat = this.lat.addPoint(ll.latRadians);
        var newLng = this.lng.addPoint(ll.lngRadians);
        return new S2LatLngRect(newLat, newLng);
    };
    /**
     * Return a rectangle that contains all points whose latitude distance from
     * this rectangle is at most margin.lat, and whose longitude distance from
     * this rectangle is at most margin.lng. In particular, latitudes are
     * clamped while longitudes are wrapped. Note that any expansion of an empty
     * interval remains empty, and both components of the given margin must be
     * non-negative.
     *
     * NOTE: If you are trying to grow a rectangle by a certain *distance* on the
     * sphere (e.g. 5km), use the ConvolveWithCap() method instead.
     */
    S2LatLngRect.prototype.expanded = function (margin) {
        // assert (margin.latRadians >= 0 && margin.lngRadians >= 0);
        if (this.isEmpty()) {
            return this;
        }
        return new S2LatLngRect(this.lat
            .expanded(margin.latRadians)
            .intersection(S2LatLngRect.fullLat()), this.lng.expanded(margin.lngRadians));
    };
    /**
     * Return the smallest rectangle containing the union of this rectangle and
     * the given rectangle.
     */
    S2LatLngRect.prototype.union = function (other) {
        return new S2LatLngRect(this.lat.union(other.lat), this.lng.union(other.lng));
    };
    /**
     * Return the smallest rectangle containing the intersection of this rectangle
     * and the given rectangle. Note that the region of intersection may consist
     * of two disjoint rectangles, in which case a single rectangle spanning both
     * of them is returned.
     */
    S2LatLngRect.prototype.intersection = function (other) {
        var intersectLat = this.lat.intersection(other.lat);
        var intersectLng = this.lng.intersection(other.lng);
        if (intersectLat.isEmpty() || intersectLng.isEmpty()) {
            // The lat/lng ranges must either be both empty or both non-empty.
            return S2LatLngRect.empty();
        }
        return new S2LatLngRect(intersectLat, intersectLng);
    };
    //
    // /**
    //  * Return a rectangle that contains the convolution of this rectangle with a
    //  * cap of the given angle. This expands the rectangle by a fixed distance (as
    //  * opposed to growing the rectangle in latitude-longitude space). The returned
    //  * rectangle includes all points whose minimum distance to the original
    //  * rectangle is at most the given angle.
    //  */
    // public S2LatLngRect convolveWithCap(/*S1Angle*/ angle) {
    //   // The most straightforward approach is to build a cap centered on each
    //   // vertex and take the union of all the bounding rectangles (including the
    //   // original rectangle; this is necessary for very large rectangles).
    //
    //   // Optimization: convert the angle to a height exactly once.
    //   S2Cap cap = S2Cap.fromAxisAngle(new S2Point(1, 0, 0), angle);
    //
    //   S2LatLngRect r = this;
    //   for (int k = 0; k < 4; ++k) {
    //     S2Cap vertexCap = S2Cap.fromAxisHeight(getVertex(k).toPoint(), cap
    //         .height());
    //     r = r.union(vertexCap.getRectBound());
    //   }
    //   return r;
    // }
    /** Return the surface area of this rectangle on the unit sphere. */
    S2LatLngRect.prototype.area = function () {
        if (this.isEmpty()) {
            return S2.toDecimal(0);
        }
        // This is the size difference of the two spherical caps, multiplied by
        // the longitude ratio.
        //TODO: check if this.lat.hi & this.lat.lo is radians. 
        return this.lng.getLength().times(decimal.Decimal.sin(this.lat.hi).minus(decimal.Decimal.sin(this.lat.lo)).abs());
    };
    /** Return true if two rectangles contains the same set of points. */
    S2LatLngRect.prototype.equals = function (that) {
        if (!(that instanceof S2LatLngRect)) {
            return false;
        }
        return this.lat.equals(that.lat) && this.lng.equals(that.lng);
    };
    /**
     * Return true if the latitude and longitude intervals of the two rectangles
     * are the same up to the given tolerance (see r1interval.h and s1interval.h
     * for details).
     */
    S2LatLngRect.prototype.approxEquals = function (other, maxError) {
        if (maxError === void 0) { maxError = 1e-15; }
        return (this.lat.approxEquals(other.lat, maxError) && this.lng.approxEquals(other.lng, maxError));
    };
    // //////////////////////////////////////////////////////////////////////
    // S2Region interface (see {@code S2Region} for details):
    S2LatLngRect.prototype.clone = function () {
        return new S2LatLngRect(this.lat, this.lng);
    };
    S2LatLngRect.prototype.getCapBound = function () {
        // We consider two possible bounding caps, one whose axis passes
        // through the center of the lat-long rectangle and one whose axis
        // is the north or south pole. We return the smaller of the two caps.
        if (this.isEmpty()) {
            return S2Cap.empty();
        }
        var poleZ, poleAngle;
        if (this.lat.lo.plus(this.lat.hi).lt(0)) {
            // South pole axis yields smaller cap.
            poleZ = -1;
            poleAngle = this.lat.hi.plus(S2.M_PI_2);
        }
        else {
            poleZ = 1;
            poleAngle = this.lat.lo.neg().plus(S2.M_PI_2);
        }
        var poleCap = S2Cap.fromAxisAngle(new S2Point(0, 0, poleZ), new S1Angle(poleAngle));
        // For bounding rectangles that span 180 degrees or less in longitude, the
        // maximum cap size is achieved at one of the rectangle vertices. For
        // rectangles that are larger than 180 degrees, we punt and always return a
        // bounding cap centered at one of the two poles.
        var lngSpan = this.lng.hi.minus(this.lng.lo);
        if (S2.IEEEremainder(lngSpan, 2 * S2.M_PI).gte(0)) {
            if (lngSpan.lt(2 * S2.M_PI)) {
                var midCap = S2Cap.fromAxisAngle(this.getCenter().toPoint(), new S1Angle(0));
                for (var k = 0; k < 4; ++k) {
                    midCap = midCap.addPoint(this.getVertex(k).toPoint());
                }
                if (midCap.height.lt(poleCap.height)) {
                    return midCap;
                }
            }
        }
        return poleCap;
    };
    S2LatLngRect.prototype.getRectBound = function () {
        return this;
    };
    S2LatLngRect.prototype.containsC = function (cell) {
        // A latitude-longitude rectangle contains a cell if and only if it contains
        // the cell's bounding rectangle. (This is an exact test.)
        return this.containsLLR(cell.getRectBound());
    };
    /**
     * This test is cheap but is NOT exact. Use Intersects() if you want a more
     * accurate and more expensive test. Note that when this method is used by an
     * S2RegionCoverer, the accuracy isn't all that important since if a cell may
     * intersect the region then it is subdivided, and the accuracy of this method
     * goes up as the cells get smaller.
     */
    S2LatLngRect.prototype.mayIntersectC = function (cell) {
        // This test is cheap but is NOT exact (see s2latlngrect.h).
        return this.intersectsLLR(cell.getRectBound());
    };
    /** The point 'p' does not need to be normalized. */
    S2LatLngRect.prototype.containsP = function (p) {
        return this.containsLL(S2LatLng.fromPoint(p));
    };
    /**
     * Return true if the edge AB intersects the given edge of constant longitude.
     */
    S2LatLngRect.intersectsLngEdge = function (a, b, lat, lng) {
        // Return true if the segment AB intersects the given edge of constant
        // longitude. The nice thing about edges of constant longitude is that
        // they are straight lines on the sphere (geodesics).
        return S2.simpleCrossing(a, b, new S2LatLng(lat.lo, lng)
            .toPoint(), new S2LatLng(lat.hi, lng).toPoint());
    };
    /**
     * Return true if the edge AB intersects the given edge of constant latitude.
     */
    S2LatLngRect.intersectsLatEdge = function (a, b, lat, lng) {
        // Return true if the segment AB intersects the given edge of constant
        // latitude. Unfortunately, lines of constant latitude are curves on
        // the sphere. They can intersect a straight edge in 0, 1, or 2 points.
        // assert (S2.isUnitLength(a) && S2.isUnitLength(b));
        // First, compute the normal to the plane AB that points vaguely north.
        var z = S2Point.normalize(S2.robustCrossProd(a, b));
        if (z.z.lt(0)) {
            z = S2Point.neg(z);
        }
        // Extend this to an orthonormal frame (x,y,z) where x is the direction
        // where the great circle through AB achieves its maximium latitude.
        var y = S2Point.normalize(S2.robustCrossProd(z, new S2Point(0, 0, 1)));
        var x = S2Point.crossProd(y, z);
        // assert (S2.isUnitLength(x) && x.z >= 0);
        // Compute the angle "theta" from the x-axis (in the x-y plane defined
        // above) where the great circle intersects the given line of latitude.
        var sinLat = decimal.Decimal.sin(lat);
        if (sinLat.abs().gte(x.z)) {
            return false; // The great circle does not reach the given latitude.
        }
        // assert (x.z > 0);
        var cosTheta = sinLat.dividedBy(x.z);
        var sinTheta = cosTheta.pow(2).neg().plus(1).sqrt(); // Math.sqrt(1 - cosTheta * cosTheta);
        var theta = decimal.Decimal.atan2(sinTheta, cosTheta);
        // Math.atan2(sinTheta, cosTheta);
        // The candidate intersection points are located +/- theta in the x-y
        // plane. For an intersection to be valid, we need to check that the
        // intersection point is contained in the interior of the edge AB and
        // also that it is contained within the given longitude interval "lng".
        // Compute the range of theta values spanned by the edge AB.
        var abTheta = S1Interval.fromPointPair(decimal.Decimal.atan2(a.dotProd(y), a.dotProd(x)), decimal.Decimal.atan2(b.dotProd(y), b.dotProd(x)));
        if (abTheta.contains(theta)) {
            // Check if the intersection point is also in the given "lng" interval.
            var isect = S2Point.add(S2Point.mul(x, cosTheta), S2Point.mul(y, sinTheta));
            if (lng.contains(decimal.Decimal.atan2(isect.y, isect.x))) {
                return true;
            }
        }
        if (abTheta.contains(theta.neg())) {
            // Check if the intersection point is also in the given "lng" interval.
            var intersection = S2Point.sub(S2Point.mul(x, cosTheta), S2Point.mul(y, sinTheta));
            if (lng.contains(decimal.Decimal.atan2(intersection.y, intersection.x))) {
                return true;
            }
        }
        return false;
    };
    S2LatLngRect.prototype.allVertex = function () {
        return [
            this.getVertex(0),
            this.getVertex(1),
            this.getVertex(2),
            this.getVertex(3)
        ];
    };
    S2LatLngRect.prototype.toGEOJSON = function () {
        return {
            type: 'Feature',
            geometry: {
                type: 'Polygon',
                coordinates: [this.allVertex().concat(this.getVertex(0)).map(function (v) { return [parseFloat(v.lngDegrees.toFixed(5)), parseFloat(v.latDegrees.toFixed(5))]; })],
            },
            properties: {}
        };
    };
    S2LatLngRect.prototype.toString = function () {
        return "[Lo=" + this.lo().toString() + ", Hi=" + this.hi().toString() + "]";
    };
    return S2LatLngRect;
}());

/*
 * Copyright 2005 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/**
 * This class represents a spherical cap, i.e. a portion of a sphere cut off by
 * a plane. The cap is defined by its axis and height. This representation has
 * good numerical accuracy for very small caps (unlike the (axis,
 * min-distance-from-origin) representation), and is also efficient for
 * containment tests (unlike the (axis, angle) representation).
 *
 * Here are some useful relationships between the cap height (h), the cap
 * opening angle (theta), the maximum chord length from the cap's center (d),
 * and the radius of cap's base (a). All formulas assume a unit radius.
 *
 * h = 1 - cos(theta) = 2 sin^2(theta/2) d^2 = 2 h = a^2 + h^2
 *
 */
var S2Cap = /** @class */ (function () {
    /**
     * Create a cap given its axis and the cap height, i.e. the maximum projected
     * distance along the cap axis from the cap center. 'axis' should be a
     * unit-length vector.
     */
    function S2Cap(axis, _height) {
        this.axis = axis;
        this.height = S2.toDecimal(_height);
        // assert (isValid());
    }
    /**
     * Create a cap given its axis and the cap opening angle, i.e. maximum angle
     * between the axis and a point on the cap. 'axis' should be a unit-length
     * vector, and 'angle' should be between 0 and 180 degrees.
     */
    S2Cap.fromAxisAngle = function (axis, angle) {
        // The height of the cap can be computed as 1-cos(angle), but this isn't
        // very accurate for angles close to zero (where cos(angle) is almost 1).
        // Computing it as 2*(sin(angle/2)**2) gives much better precision.
        // assert (S2.isUnitLength(axis));
        var d = angle.radians.times(0.5).sin();
        // ecimal.sin(0.5 * angle.radians.times(0.5));
        return new S2Cap(axis, d.pow(2).times(2));
    };
    /**
     * Create a cap given its axis and its area in steradians. 'axis' should be a
     * unit-length vector, and 'area' should be between 0 and 4 * M_PI.
     */
    S2Cap.fromAxisArea = function (axis, _area) {
        var area = S2.toDecimal(_area);
        // assert (S2.isUnitLength(axis));
        return new S2Cap(axis, area.dividedBy(S2.toDecimal(2).times(S2.M_PI)));
    };
    /** Return an empty cap, i.e. a cap that contains no points. */
    S2Cap.empty = function () {
        return new S2Cap(new S2Point(1, 0, 0), -1);
    };
    /** Return a full cap, i.e. a cap that contains all points. */
    S2Cap.full = function () {
        return new S2Cap(new S2Point(1, 0, 0), 2);
    };
    S2Cap.prototype.getCapBound = function () {
        return this;
    };
    S2Cap.prototype.area = function () {
        return decimal.Decimal.max(0, this.height)
            .times(S2.M_PI)
            .times(2);
        // return 2 * S2.M_PI * Math.max(0.0, this.height);
    };
    /**
     * Return the cap opening angle in radians, or a negative number for empty
     * caps.
     */
    S2Cap.prototype.angle = function () {
        // This could also be computed as acos(1 - height_), but the following
        // formula is much more accurate when the cap height is small. It
        // follows from the relationship h = 1 - cos(theta) = 2 sin^2(theta/2).
        if (this.isEmpty()) {
            return new S1Angle(-1);
        }
        return new S1Angle(decimal.Decimal.asin(this.height.times(0.5).sqrt())
            .times(2));
    };
    /**
     * We allow negative heights (to represent empty caps) but not heights greater
     * than 2.
     */
    S2Cap.prototype.isValid = function () {
        return S2.isUnitLength(this.axis) && this.height.lte(2);
    };
    /** Return true if the cap is empty, i.e. it contains no points. */
    S2Cap.prototype.isEmpty = function () {
        return this.height.lt(0);
    };
    /** Return true if the cap is full, i.e. it contains all points. */
    S2Cap.prototype.isFull = function () {
        return this.height.gte(2);
    };
    /**
     * Return the complement of the interior of the cap. A cap and its complement
     * have the same boundary but do not share any interior points. The complement
     * operator is not a bijection, since the complement of a singleton cap
     * (containing a single point) is the same as the complement of an empty cap.
     */
    S2Cap.prototype.complement = function () {
        // The complement of a full cap is an empty cap, not a singleton.
        // Also make sure that the complement of an empty cap has height 2.
        var cHeight = this.isFull() ? -1 : decimal.Decimal.max(this.height, 0).neg().plus(2);
        return new S2Cap(S2Point.neg(this.axis), cHeight);
    };
    /**
     * Return true if and only if this cap contains the given other cap (in a set
     * containment sense, e.g. every cap contains the empty cap).
     */
    S2Cap.prototype.containsCap = function (other) {
        if (this.isFull() || other.isEmpty()) {
            return true;
        }
        return this.angle().radians.gte(this.axis.angle(other.axis).plus(other.angle().radians));
    };
    /**
     * Return true if and only if the interior of this cap intersects the given
     * other cap. (This relationship is not symmetric, since only the interior of
     * this cap is used.)
     */
    S2Cap.prototype.interiorIntersects = function (other) {
        // Interior(X) intersects Y if and only if Complement(Interior(X))
        // does not contain Y.
        return !this.complement().containsCap(other);
    };
    /**
     * Return true if and only if the given point is contained in the interior of
     * the region (i.e. the region excluding its boundary). 'p' should be a
     * unit-length vector.
     */
    S2Cap.prototype.interiorContains = function (p) {
        // assert (S2.isUnitLength(p));
        return this.isFull() || S2Point.sub(this.axis, p).norm2().lt(this.height.times(2));
    };
    /**
     * Increase the cap height if necessary to include the given point. If the cap
     * is empty the axis is set to the given point, but otherwise it is left
     * unchanged. 'p' should be a unit-length vector.
     */
    S2Cap.prototype.addPoint = function (p) {
        // Compute the squared chord length, then convert it into a height.
        // assert (S2.isUnitLength(p));
        if (this.isEmpty()) {
            return new S2Cap(p, 0);
        }
        else {
            // To make sure that the resulting cap actually includes this point,
            // we need to round up the distance calculation. That is, after
            // calling cap.AddPoint(p), cap.Contains(p) should be true.
            var dist2 = S2Point.sub(this.axis, p).norm2();
            var newHeight = decimal.Decimal.max(this.height, S2Cap.ROUND_UP.times(0.5).times(dist2));
            return new S2Cap(this.axis, newHeight);
        }
    };
    // Increase the cap height if necessary to include "other". If the current
    // cap is empty it is set to the given other cap.
    S2Cap.prototype.addCap = function (other) {
        if (this.isEmpty()) {
            return new S2Cap(other.axis, other.height);
        }
        else {
            // See comments for FromAxisAngle() and AddPoint(). This could be
            // optimized by doing the calculation in terms of cap heights rather
            // than cap opening angles.
            var angle = this.axis.angle(other.axis).plus(other.angle().radians);
            if (angle.gte(S2.M_PI)) {
                return new S2Cap(this.axis, 2); //Full cap
            }
            else {
                var d = angle.times(0.5).sin();
                var newHeight = decimal.Decimal.max(this.height, S2Cap.ROUND_UP.times(2).times(d.pow(2)));
                return new S2Cap(this.axis, newHeight);
            }
        }
    };
    // //////////////////////////////////////////////////////////////////////
    // S2Region interface (see {@code S2Region} for details):
    S2Cap.prototype.getRectBound = function () {
        if (this.isEmpty()) {
            return S2LatLngRect.empty();
        }
        // Convert the axis to a (lat,lng) pair, and compute the cap angle.
        var axisLatLng = S2LatLng.fromPoint(this.axis);
        var capAngle = this.angle().radians;
        var allLongitudes = false;
        var lat = Array(2);
        var lng = Array(2);
        lng[0] = S2.toDecimal(-S2.M_PI);
        lng[1] = S2.toDecimal(S2.M_PI);
        // Check whether cap includes the south pole.
        lat[0] = axisLatLng.latRadians.minus(capAngle);
        if (lat[0].lte(-S2.M_PI_2)) {
            lat[0] = S2.toDecimal(-S2.M_PI_2);
            allLongitudes = true;
        }
        // Check whether cap includes the north pole.
        lat[1] = axisLatLng.latRadians.plus(capAngle);
        if (lat[1].gte(S2.M_PI_2)) {
            lat[1] = S2.toDecimal(S2.M_PI_2);
            allLongitudes = true;
        }
        if (!allLongitudes) {
            // Compute the range of longitudes covered by the cap. We use the law
            // of sines for spherical triangles. Consider the triangle ABC where
            // A is the north pole, B is the center of the cap, and C is the point
            // of tangency between the cap boundary and a line of longitude. Then
            // C is a right angle, and letting a,b,c denote the sides opposite A,B,C,
            // we have sin(a)/sin(A) = sin(c)/sin(C), or sin(A) = sin(a)/sin(c).
            // Here "a" is the cap angle, and "c" is the colatitude (90 degrees
            // minus the latitude). This formula also works for negative latitudes.
            //
            // The formula for sin(a) follows from the relationship h = 1 - cos(a).
            // double sinA = Math.sqrt(this.height * (2 - this.height));
            // double sinC = Math.cos(axisLatLng.lat().radians());
            var sinA = this.height.times(this.height.neg().plus(2)).sqrt();
            var sinC = axisLatLng.latRadians.cos();
            if (sinA.lte(sinC)) {
                var angleA = decimal.Decimal.asin(sinA.dividedBy(sinC));
                lng[0] = S2.IEEEremainder(axisLatLng.lngRadians.minus(angleA), 2 * S2.M_PI);
                lng[1] = S2.IEEEremainder(axisLatLng.lngRadians.plus(angleA), 2 * S2.M_PI);
            }
        }
        return new S2LatLngRect(new R1Interval(lat[0], lat[1]), new S1Interval(lng[0], lng[1]));
    };
    S2Cap.prototype.containsC = function (cell) {
        // If the cap does not contain all cell vertices, return false.
        // We check the vertices before taking the Complement() because we can't
        // accurately represent the complement of a very small cap (a height
        // of 2-epsilon is rounded off to 2).
        var vertices = new Array(4);
        for (var k = 0; k < 4; ++k) {
            vertices[k] = cell.getVertex(k);
            if (!this.contains(vertices[k])) {
                return false;
            }
        }
        // Otherwise, return true if the complement of the cap does not intersect
        // the cell. (This test is slightly conservative, because technically we
        // want Complement().InteriorIntersects() here.)
        return !this.complement().intersects(cell, vertices);
    };
    // public mayIntersectC(cell:S2Cell):boolean {
    //   const toRet = this._mayIntersectC(cell);
    //   console.log("intersects? ",toRet, cell.id.pos().toString(16), cell.level);
    //   return toRet;
    // }
    S2Cap.prototype.mayIntersectC = function (cell) {
        // If the cap contains any cell vertex, return true.
        var vertices = new Array(4);
        for (var k = 0; k < 4; ++k) {
            vertices[k] = cell.getVertex(k);
            if (this.contains(vertices[k])) {
                return true;
            }
        }
        return this.intersects(cell, vertices);
    };
    /**
     * Return true if the cap intersects 'cell', given that the cap vertices have
     * alrady been checked.
     */
    S2Cap.prototype.intersects = function (cell, vertices) {
        // Return true if this cap intersects any point of 'cell' excluding its
        // vertices (which are assumed to already have been checked).
        // If the cap is a hemisphere or larger, the cell and the complement of the
        // cap are both convex. Therefore since no vertex of the cell is contained,
        // no other interior point of the cell is contained either.
        if (this.height.gte(1)) {
            return false;
        }
        // We need to check for empty caps due to the axis check just below.
        if (this.isEmpty()) {
            return false;
        }
        // Optimization: return true if the cell contains the cap axis. (This
        // allows half of the edge checks below to be skipped.)
        if (cell.contains(this.axis)) {
            return true;
        }
        // At this point we know that the cell does not contain the cap axis,
        // and the cap does not contain any cell vertex. The only way that they
        // can intersect is if the cap intersects the interior of some edge.
        var sin2Angle = this.height.times(this.height.neg().plus(2)); // sin^2(capAngle)
        // if (cell.id.pos().toString(16) === '77c040000000000') {
        //   console.log("DIOCAN");
        // }
        for (var k = 0; k < 4; ++k) {
            var edge = cell.getEdgeRaw(k);
            var dot = this.axis.dotProd(edge);
            if (dot.gt(0)) {
                // The axis is in the interior half-space defined by the edge. We don't
                // need to consider these edges, since if the cap intersects this edge
                // then it also intersects the edge on the opposite side of the cell
                // (because we know the axis is not contained with the cell).
                continue;
            }
            // The Norm2() factor is necessary because "edge" is not normalized.
            if (dot.pow(2).gt(sin2Angle.times(edge.norm2()))) {
                // if (cell.id.pos().toString(16) === '77c040000000000') {
                //   console.log("DIOCaAN", k, dot.toString(), sin2Angle.toString(), sin2Angle.times(edge.norm2()).toString());
                // }
                return false; // Entire cap is on the exterior side of this edge.
            }
            // Otherwise, the great circle containing this edge intersects
            // the interior of the cap. We just need to check whether the point
            // of closest approach occurs between the two edge endpoints.
            var dir = S2Point.crossProd(edge, this.axis);
            if (dir.dotProd(vertices[k]).lt(0)
                && dir.dotProd(vertices[(k + 1) & 3]).gt(0)) {
                return true;
            }
        }
        return false;
    };
    S2Cap.prototype.contains = function (p) {
        // The point 'p' should be a unit-length vector.
        // assert (S2.isUnitLength(p));
        return S2Point.sub(this.axis, p).norm2().lte(this.height.times(2));
    };
    //
    // /** Return true if two caps are identical. */
    // public equals(that:Object ):boolean  {
    //
    //   if (!(that instanceof S2Cap)) {
    //     return false;
    //   }
    //
    //   S2Cap other = (S2Cap) that;
    //   return (this.axis.equals(other.axis) && this.height == other.height)
    //       || (isEmpty() && other.isEmpty()) || (isFull() && other.isFull());
    //
    // }
    //
    // @Override
    // public int hashCode() {
    //   if (isFull()) {
    //     return 17;
    //   } else if (isEmpty()) {
    //     return 37;
    //   }
    //   int result = 17;
    //   result = 37 * result + this.axis.hashCode();
    //   long heightBits = Double.doubleToLongBits(this.height);
    //   result = 37 * result + (int) ((heightBits >>> 32) ^ heightBits);
    //   return result;
    // }
    // /////////////////////////////////////////////////////////////////////
    // The following static methods are convenience functions for assertions
    // and testing purposes only.
    /**
     * Return true if the cap axis and height differ by at most "max_error" from
     * the given cap "other".
     */
    S2Cap.prototype.approxEquals = function (other, maxError) {
        if (maxError === void 0) { maxError = 1e-14; }
        return (this.axis.aequal(other.axis, maxError) && this.height.minus(other.height).lte(maxError))
            || (this.isEmpty() && other.height.lte(maxError))
            || (other.isEmpty() && this.height.lte(maxError))
            || (this.isFull() && other.height.gte(2 - maxError))
            || (other.isFull() && this.height.gte(2 - maxError));
    };
    S2Cap.prototype.toString = function () {
        return "[Point = " + this.axis.toString() + " Height = " + this.height.toString() + "]";
    };
    S2Cap.prototype.toGEOJSON = function () {
        return this.getRectBound().toGEOJSON();
    };
    /**
     * Multiply a positive number by this constant to ensure that the result of a
     * floating point operation is at least as large as the true
     * infinite-precision result.
     */
    S2Cap.ROUND_UP = S2.toDecimal(1).dividedBy(new Long(1).shiftLeft(52).toString()).plus(1);
    return S2Cap;
}());

var MutableInteger = /** @class */ (function () {
    function MutableInteger(val) {
        this.val = val;
    }
    return MutableInteger;
}());

/*
 * Copyright 2005 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
var parseHex = function parseHex(str) {
    return Long.fromString(str, false, 16);
};
/**
 * An S2CellId is a 64-bit unsigned integer that uniquely identifies a cell in
 * the S2 cell decomposition. It has the following format:
 *
 * <pre>
 * id = [face][face_pos]
 * </pre>
 *
 * face: a 3-bit number (range 0..5) encoding the cube face.
 *
 * face_pos: a 61-bit number encoding the position of the center of this cell
 * along the Hilbert curve over this face (see the Wiki pages for details).
 *
 * Sequentially increasing cell ids follow a continuous space-filling curve over
 * the entire sphere. They have the following properties:
 *  - The id of a cell at level k consists of a 3-bit face number followed by k
 * bit pairs that recursively select one of the four children of each cell. The
 * next bit is always 1, and all other bits are 0. Therefore, the level of a
 * cell is determined by the position of its lowest-numbered bit that is turned
 * on (for a cell at level k, this position is 2 * (MAX_LEVEL - k).)
 *  - The id of a parent cell is at the midpoint of the range of ids spanned by
 * its children (or by its descendants at any level).
 *
 * Leaf cells are often used to represent points on the unit sphere, and this
 * class provides methods for converting directly between these two
 * representations. For cells that represent 2D regions rather than discrete
 * point, it is better to use the S2Cell class.
 *
 *
 */
var S2CellId = /** @class */ (function () {
    function S2CellId(id) {
        if (typeof (id) === 'string') {
            this.id = Long.fromString(id);
        }
        else {
            this.id = id;
        }
    }
    Object.defineProperty(S2CellId.prototype, "face", {
        /** Which cube face this cell belongs to, in the range 0..5. */
        get: function () {
            return this.id.shiftRightUnsigned(S2CellId.POS_BITS).toInt();
        },
        enumerable: true,
        configurable: true
    });
    /** Return the lowest-numbered bit that is on for cells at the given level. */
    S2CellId.prototype.lowestOnBit = function () {
        return this.id.and(this.id.negate());
    };
    /** The default constructor returns an invalid cell id. */
    S2CellId.none = function () {
        return new S2CellId(new Long(0));
    };
    /**
     * Returns an invalid cell id guaranteed to be larger than any valid cell id.
     * Useful for creating indexes.
     */
    S2CellId.sentinel = function () {
        return new S2CellId(S2CellId.MAX_UNSIGNED); // -1
    };
    S2CellId.prototype.getBits1 = function (i, j, k, bits) {
        var nbits = (k == 7) ? (S2CellId.MAX_LEVEL - 7 * S2CellId.LOOKUP_BITS) : S2CellId.LOOKUP_BITS;
        bits += (this.id
            .shiftRightUnsigned((k * 2 * S2CellId.LOOKUP_BITS + 1))
            .getLowBitsUnsigned()
            & ((1 << (2 * nbits)) - 1)) << 2;
        /*
         * System.out.println("id is: " + id_); System.out.println("bits is " +
         * bits); System.out.println("lookup_ij[bits] is " + lookup_ij[bits]);
         */
        bits = S2CellId.LOOKUP_IJ[bits];
        i.val = i.val + ((bits >> (S2CellId.LOOKUP_BITS + 2)) << (k * S2CellId.LOOKUP_BITS));
        // i.setValue(i.intValue() + ((bits >> (LOOKUP_BITS + 2)) << (k * LOOKUP_BITS)));
        /*
         * System.out.println("left is " + ((bits >> 2) & ((1 << kLookupBits) -
         * 1))); System.out.println("right is " + (k * kLookupBits));
         * System.out.println("j is: " + j.intValue()); System.out.println("addition
         * is: " + ((((bits >> 2) & ((1 << kLookupBits) - 1))) << (k *
         * kLookupBits)));
         */
        j.val = j.val + ((((bits >> 2) & ((1 << S2CellId.LOOKUP_BITS) - 1))) << (k * S2CellId.LOOKUP_BITS));
        bits &= (S2.SWAP_MASK | S2.INVERT_MASK);
        return bits;
    };
    /**
     * Convert (face, si, ti) coordinates (see s2.h) to a direction vector (not
     * necessarily unit length).
     */
    S2CellId.prototype.faceSiTiToXYZ = function (face, si, ti) {
        // console.log('faceSiTiToXYZ', si, ti);
        var kScale = S2.toDecimal(1).dividedBy(S2CellId.MAX_SIZE);
        var uvVector = R2Vector.fromSTVector(new R2Vector(kScale.times(si), kScale.times(ti)));
        // console.log(uvVector.toString(), uvVector.x.toString());
        return uvVector.toPoint(face);
    };
    S2CellId.lowestOnBitForLevel = function (level) {
        return new Long(1).shiftLeft(2 * (S2CellId.MAX_LEVEL - level));
    };
    /**
     * Return the (face, i, j) coordinates for the leaf cell corresponding to this
     * cell id. Since cells are represented by the Hilbert curve position at the
     * center of the cell, the returned (i,j) for non-leaf cells will be a leaf
     * cell adjacent to the cell center. If "orientation" is non-NULL, also return
     * the Hilbert curve orientation for the current cell.
     */
    S2CellId.prototype.toFaceIJOrientation = function (pi, pj, orientation) {
        // System.out.println("Entering toFaceIjorientation");
        var face = this.face;
        var bits = (face & S2.SWAP_MASK);
        // System.out.println("face = " + face + " bits = " + bits);
        // Each iteration maps 8 bits of the Hilbert curve position into
        // 4 bits of "i" and "j". The lookup table transforms a key of the
        // form "ppppppppoo" to a value of the form "iiiijjjjoo", where the
        // letters [ijpo] represents bits of "i", "j", the Hilbert curve
        // position, and the Hilbert curve orientation respectively.
        //
        // On the first iteration we need to be careful to clear out the bits
        // representing the cube face.
        for (var k = 7; k >= 0; --k) {
            bits = this.getBits1(pi, pj, k, bits);
            // System.out.println("pi = " + pi + " pj= " + pj + " bits = " + bits);
        }
        if (orientation != null) {
            // The position of a non-leaf cell at level "n" consists of a prefix of
            // 2*n bits that identifies the cell, followed by a suffix of
            // 2*(MAX_LEVEL-n)+1 bits of the form 10*. If n==MAX_LEVEL, the suffix is
            // just "1" and has no effect. Otherwise, it consists of "10", followed
            // by (MAX_LEVEL-n-1) repetitions of "00", followed by "0". The "10" has
            // no effect, while each occurrence of "00" has the effect of reversing
            // the kSwapMask bit.
            // assert (S2.POS_TO_ORIENTATION[2] == 0);
            // assert (S2.POS_TO_ORIENTATION[0] == S2.SWAP_MASK);
            if ((Long.fromString('0x1111111111111110', true, 16).and(this.lowestOnBit()).notEquals(0))) {
                bits ^= S2.SWAP_MASK;
            }
            orientation.val = bits;
        }
        return face;
    };
    /**
     * Return true if this is a leaf cell (more efficient than checking whether
     * level() == MAX_LEVEL).
     */
    S2CellId.prototype.isLeaf = function () {
        return this.id.and(1).getLowBits() != 0;
    };
    /**
     * Return the cell at the previous level or at the given level (which must be
     * less than or equal to the current level).
     */
    S2CellId.prototype.parentL = function (level) {
        // assert (isValid() && level >= 0 && level <= this.level());
        var newLsb = S2CellId.lowestOnBitForLevel(level);
        return new S2CellId(this.id.and(newLsb.negate()).or(newLsb));
        // return new S2CellId((id & -newLsb) | newLsb);
    };
    S2CellId.prototype.parent = function () {
        // assert (isValid() && level() > 0);
        var newLsb = this.lowestOnBit().shiftLeft(2);
        // return new S2CellId((id & -newLsb) | newLsb);
        return new S2CellId(this.id.and(newLsb.negate()).or(newLsb));
    };
    /**
     * Return a cell given its face (range 0..5), 61-bit Hilbert curve position
     * within that face, and level (range 0..MAX_LEVEL). The given position will
     * be modified to correspond to the Hilbert curve position at the center of
     * the returned cell. This is a static function rather than a constructor in
     * order to give names to the arguments.
     */
    S2CellId.fromFacePosLevel = function (face, pos, level) {
        // equivalent to pos | 1
        return new S2CellId(new Long(face)
            .shiftLeft(S2CellId.POS_BITS)
            .add(pos.or(1))).parentL(level);
        // return new S2CellId((((long) face) << POS_BITS) + (pos | 1)).parent(level);
    };
    // /**
    //  * Return the leaf cell containing the given point (a direction vector, not
    //  * necessarily unit length).
    //  */
    S2CellId.fromPoint = function (p) {
        var face = p.toFace();
        var uv = p.toR2Vector(face);
        var i = S2CellId.stToIJ(uv.toSt(0));
        var j = S2CellId.stToIJ(uv.toSt(1));
        return S2CellId.fromFaceIJ(face, i, j);
    };
    //
    //
    // /** Return the leaf cell containing the given S2LatLng. */
    // public static S2CellId fromLatLng(S2LatLng ll) {
    //   return fromPoint(ll.toPoint());
    // }
    S2CellId.prototype.toPoint = function () {
        return S2Point.normalize(this.toPointRaw());
    };
    /**
     * Return the direction vector corresponding to the center of the given cell.
     * The vector returned by ToPointRaw is not necessarily unit length.
     */
    S2CellId.prototype.toPointRaw = function () {
        // First we compute the discrete (i,j) coordinates of a leaf cell contained
        // within the given cell. Given that cells are represented by the Hilbert
        // curve position corresponding at their center, it turns out that the cell
        // returned by ToFaceIJOrientation is always one of two leaf cells closest
        // to the center of the cell (unless the given cell is a leaf cell itself,
        // in which case there is only one possibility).
        //
        // Given a cell of size s >= 2 (i.e. not a leaf cell), and letting (imin,
        // jmin) be the coordinates of its lower left-hand corner, the leaf cell
        // returned by ToFaceIJOrientation() is either (imin + s/2, jmin + s/2)
        // (imin + s/2 - 1, jmin + s/2 - 1). We can distinguish these two cases by
        // looking at the low bit of "i" or "j". In the first case the low bit is
        // zero, unless s == 2 (i.e. the level just above leaf cells) in which case
        // the low bit is one.
        //
        // The following calculation converts (i,j) to the (si,ti) coordinates of
        // the cell center. (We need to multiply the coordinates by a factor of 2
        // so that the center of leaf cells can be represented exactly.)
        var i = new MutableInteger(0);
        var j = new MutableInteger(0);
        var face = this.toFaceIJOrientation(i, j, null);
        // System.out.println("i= " + i.intValue() + " j = " + j.intValue());
        // let delta = isLeaf() ? 1 : (((i.intValue() ^ (((int) id) >>> 2)) & 1) != 0) ? 2 : 0;
        var delta = this.isLeaf()
            ? 1 :
            ((((new Long(i.val).getLowBits() ^ ((this.id.getLowBits()) >>> 2)) & 1) != 0)
                ? 2 : 0);
        // let delta = this.isLeaf() ? 1 : new Long(i.val).and(this.id.getLowBits() >>> 2).and(1).notEquals(1) ? 2 : 0
        // ((i.val ? (((int)id) >>> 2))  & 1  ))
        var si = new Long((i.val << 1) + delta - S2CellId.MAX_SIZE).getLowBits();
        var ti = new Long((j.val << 1) + delta - S2CellId.MAX_SIZE).getLowBits();
        return this.faceSiTiToXYZ(face, si, ti);
    };
    /** Return the S2LatLng corresponding to the center of the given cell. */
    S2CellId.prototype.toLatLng = function () {
        return S2LatLng.fromPoint(this.toPointRaw());
    };
    /** Return true if id() represents a valid cell. */
    S2CellId.prototype.isValid = function () {
        return this.face < S2CellId.NUM_FACES && ((this.lowestOnBit().and(Long.fromString('0x1555555555555555', false, 16)).notEquals(0)));
        // return this.face() < NUM_FACES && ((lowestOnBit() & (0x1555555555555555L)) != 0);
    };
    /**
     * The position of the cell center along the Hilbert curve over this face, in
     * the range 0..(2**kPosBits-1).
     */
    S2CellId.prototype.pos = function () {
        return this.id.and(S2CellId.MAX_UNSIGNED.shiftRightUnsigned(S2CellId.FACE_BITS));
        // return (id & (-1L >>> FACE_BITS));
    };
    /** Return the subdivision level of the cell (range 0..MAX_LEVEL). */
    S2CellId.prototype.level = function () {
        // Fast path for leaf cells.
        if (this.isLeaf()) {
            return S2CellId.MAX_LEVEL;
        }
        var x = this.id.getLowBits();
        var level = -1;
        if (x != 0) {
            level += 16;
        }
        else {
            x = this.id.shiftRightUnsigned(32).getLowBits();
            // (int) (id >>> 32);
        }
        // We only need to look at even-numbered bits to determine the
        // level of a valid cell id.
        x &= -x; // Get lowest bit.
        if ((x & 0x00005555) != 0) {
            level += 8;
        }
        if ((x & 0x00550055) != 0) {
            level += 4;
        }
        if ((x & 0x05050505) != 0) {
            level += 2;
        }
        if ((x & 0x11111111) != 0) {
            level += 1;
        }
        // assert (level >= 0 && level <= MAX_LEVEL);
        return level;
    };
    /**
     * Return true if this is a top-level face cell (more efficient than checking
     * whether level() == 0).
     */
    S2CellId.prototype.isFace = function () {
        return this.level() === 0;
        // return (id & (lowestOnBitForLevel(0) - 1)) == 0;
    };
    /**
     * Return the child position (0..3) of this cell's ancestor at the given
     * level, relative to its parent. The argument should be in the range
     * 1..MAX_LEVEL. For example, child_position(1) returns the position of this
     * cell's level-1 ancestor within its top-level face cell.
     */
    S2CellId.prototype.childPosition = function (level) {
        return this.id.shiftRight((2 * (S2CellId.MAX_LEVEL - level) + 1)).and(3).getLowBits();
        // return (int) (id >>> (2 * (MAX_LEVEL - level) + 1)) & 3;
    };
    // Methods that return the range of cell ids that are contained
    // within this cell (including itself). The range is *inclusive*
    // (i.e. test using >= and <=) and the return values of both
    // methods are valid leaf cell ids.
    //
    // These methods should not be used for iteration. If you want to
    // iterate through all the leaf cells, call child_begin(MAX_LEVEL) and
    // child_end(MAX_LEVEL) instead.
    //
    // It would in fact be error-prone to define a range_end() method,
    // because (range_max().id() + 1) is not always a valid cell id, and the
    // iterator would need to be tested using "<" rather that the usual "!=".
    S2CellId.prototype.rangeMin = function () {
        return new S2CellId(this.id.sub(this.lowestOnBit().sub(1)));
        // return new S2CellId(id - (lowestOnBit() - 1));
    };
    S2CellId.prototype.rangeMax = function () {
        return new S2CellId(this.id.add(this.lowestOnBit().sub(1)));
        // return new S2CellId(id + (lowestOnBit() - 1));
    };
    //
    //
    /** Return true if the given cell is contained within this one. */
    S2CellId.prototype.contains = function (other) {
        // assert (isValid() && other.isValid());
        return other.greaterOrEquals(this.rangeMin()) && other.lessOrEquals(this.rangeMax());
    };
    /** Return true if the given cell intersects this one. */
    S2CellId.prototype.intersects = function (other) {
        // assert (isValid() && other.isValid());
        return other.rangeMin().lessOrEquals(this.rangeMax())
            && other.rangeMax().greaterOrEquals(this.rangeMin());
    };
    S2CellId.prototype.childBegin = function () {
        // assert (isValid() && level() < MAX_LEVEL);
        var oldLsb = this.lowestOnBit();
        return new S2CellId(this.id.sub(oldLsb).add(oldLsb.shiftRight(2)));
        // return new S2CellId(id - oldLsb + (oldLsb >>> 2));
    };
    S2CellId.prototype.childBeginL = function (level) {
        // assert (isValid() && level >= this.level() && level <= MAX_LEVEL);
        return new S2CellId(this.id.sub(this.lowestOnBit()).add(S2CellId.lowestOnBitForLevel(level)));
        // return new S2CellId(id - lowestOnBit() + lowestOnBitForLevel(level));
    };
    S2CellId.prototype.childEnd = function () {
        // assert (isValid() && level() < MAX_LEVEL);
        var oldLsb = this.lowestOnBit();
        return new S2CellId(this.id.add(oldLsb).add(oldLsb.shiftRightUnsigned(2)));
        // return new S2CellId(id + oldLsb + (oldLsb >>> 2));
    };
    S2CellId.prototype.childEndL = function (level) {
        // assert (isValid() && level >= this.level() && level <= MAX_LEVEL);
        return new S2CellId(this.id.add(this.lowestOnBit()).add(S2CellId.lowestOnBitForLevel(level)));
        // return new S2CellId(id + lowestOnBit() + lowestOnBitForLevel(level));
    };
    //
    // Iterator-style methods for traversing the immediate children of a cell or
    // all of the children at a given level (greater than or equal to the current
    // level). Note that the end value is exclusive, just like standard STL
    // iterators, and may not even be a valid cell id. You should iterate using
    // code like this:
    //
    // for(S2CellId c = id.childBegin(); !c.equals(id.childEnd()); c = c.next())
    // ...
    //
    // The convention for advancing the iterator is "c = c.next()", so be sure
    // to use 'equals()' in the loop guard, or compare 64-bit cell id's,
    // rather than "c != id.childEnd()".
    /**
     * Return the next cell at the same level along the Hilbert curve. Works
     * correctly when advancing from one face to the next, but does *not* wrap
     * around from the last face to the first or vice versa.
     */
    S2CellId.prototype.next = function () {
        return new S2CellId(this.id.add(this.lowestOnBit().shiftLeft(1)));
        // return new S2CellId(id + (lowestOnBit() << 1));
    };
    /**
     * Return the previous cell at the same level along the Hilbert curve. Works
     * correctly when advancing from one face to the next, but does *not* wrap
     * around from the last face to the first or vice versa.
     */
    S2CellId.prototype.prev = function () {
        return new S2CellId(this.id.sub(this.lowestOnBit().shiftLeft(1)));
        // return new S2CellId(id - (lowestOnBit() << 1));
    };
    /**
     * Like next(), but wraps around from the last face to the first and vice
     * versa. Should *not* be used for iteration in conjunction with
     * child_begin(), child_end(), Begin(), or End().
     */
    S2CellId.prototype.nextWrap = function () {
        var n = this.next();
        if (S2CellId.unsignedLongLessThan(n.id, S2CellId.WRAP_OFFSET)) {
            return n;
        }
        return new S2CellId(n.id.sub(S2CellId.WRAP_OFFSET));
        // return new S2CellId(n.id - WRAP_OFFSET);
    };
    /**
     * Like prev(), but wraps around from the last face to the first and vice
     * versa. Should *not* be used for iteration in conjunction with
     * child_begin(), child_end(), Begin(), or End().
     */
    S2CellId.prototype.prevWrap = function () {
        var p = this.prev();
        if (p.id.lessThan(S2CellId.WRAP_OFFSET)) {
            return p;
        }
        return new S2CellId(p.id.add(S2CellId.WRAP_OFFSET));
    };
    S2CellId.begin = function (level) {
        return S2CellId.fromFacePosLevel(0, new Long(0), 0).childBeginL(level);
    };
    S2CellId.end = function (level) {
        return S2CellId.fromFacePosLevel(5, new Long(0), 0).childEndL(level);
    };
    /**
     * Decodes the cell id from a compact text string suitable for display or
     * indexing. Cells at lower levels (i.e. larger cells) are encoded into
     * fewer characters. The maximum token length is 16.
     *
     * @param token the token to decode
     * @return the S2CellId for that token
     * @throws NumberFormatException if the token is not formatted correctly
     */
    S2CellId.fromToken = function (token) {
        if (token == null) {
            throw new Error("Null string in S2CellId.fromToken");
        }
        if (token.length == 0) {
            throw new Error("Empty string in S2CellId.fromToken");
        }
        if (token.length > 16 || "X" == token) {
            return S2CellId.none();
        }
        var value = new Long(0);
        for (var pos = 0; pos < 16; pos++) {
            var digit = new Long(0);
            if (pos < token.length) {
                digit = Long.fromString(token[pos], true, 16);
                if (digit.equals(-1)) {
                    throw new Error(token);
                }
                if (S2CellId.overflowInParse(value, digit.toNumber())) {
                    throw new Error("Too large for unsigned long: " + token);
                }
            }
            value = value.mul(16).add(digit);
            // (value * 16) + digit;
        }
        return new S2CellId(value);
    };
    /**
     * Encodes the cell id to compact text strings suitable for display or indexing.
     * Cells at lower levels (i.e. larger cells) are encoded into fewer characters.
     * The maximum token length is 16.
     *
     * Simple implementation: convert the id to hex and strip trailing zeros. We
     * could use base-32 or base-64, but assuming the cells used for indexing
     * regions are at least 100 meters across (level 16 or less), the savings
     * would be at most 3 bytes (9 bytes hex vs. 6 bytes base-64).
     *
     * @return the encoded cell id
     */
    S2CellId.prototype.toToken = function () {
        if (this.id.equals(0)) {
            return "X";
        }
        var hex = this.id.toUnsigned().toString(16);
        // Long.toHexString(id).toLowerCase(Locale.ENGLISH);
        var sb = '';
        for (var i = hex.length; i < 16; i++) {
            sb += '0';
            // sb.append('0');
        }
        sb += hex;
        // sb.append(hex);
        for (var len = 16; len > 0; len--) {
            if (sb[len - 1] != '0') {
                return sb.substring(0, len);
            }
        }
        throw new Error("Shouldn't make it here");
    };
    /**
     * Returns true if (current * radix) + digit is a number too large to be
     * represented by an unsigned long.  This is useful for detecting overflow
     * while parsing a string representation of a number.
     * Does not verify whether supplied radix is valid, passing an invalid radix
     * will give undefined results or an ArrayIndexOutOfBoundsException.
     */
    S2CellId.overflowInParse = function (current, digit, radix) {
        if (radix === void 0) { radix = 10; }
        if (current.greaterThanOrEqual(0)) {
            if (current.lessThan(S2CellId.maxValueDivs[radix])) {
                return false;
            }
            if (current.greaterThan(S2CellId.maxValueDivs[radix])) {
                return true;
            }
            // current == maxValueDivs[radix]
            return (digit > S2CellId.maxValueMods[radix]);
        }
        // current < 0: high bit is set
        return true;
    };
    /**
     * Return the four cells that are adjacent across the cell's four edges.
     * Neighbors are returned in the order defined by S2Cell::GetEdge. All
     * neighbors are guaranteed to be distinct.
     */
    S2CellId.prototype.getEdgeNeighbors = function () {
        var i = new MutableInteger(0);
        var j = new MutableInteger(0);
        var level = this.level();
        var size = 1 << (S2CellId.MAX_LEVEL - level);
        var face = this.toFaceIJOrientation(i, j, null);
        var neighbors = [];
        // Edges 0, 1, 2, 3 are in the S, E, N, W directions.
        neighbors.push(S2CellId.fromFaceIJSame(face, i.val, j.val - size, j.val - size >= 0).parentL(level));
        neighbors.push(S2CellId.fromFaceIJSame(face, i.val + size, j.val, i.val + size < S2CellId.MAX_SIZE).parentL(level));
        neighbors.push(S2CellId.fromFaceIJSame(face, i.val, j.val + size, j.val + size < S2CellId.MAX_SIZE).parentL(level));
        neighbors.push(S2CellId.fromFaceIJSame(face, i.val - size, j.val, i.val - size >= 0).parentL(level));
        // neighbors[0] = fromFaceIJSame(face, i.intValue(), j.intValue() - size,
        //     j.intValue() - size >= 0).parent(level);
        // neighbors[1] = fromFaceIJSame(face, i.intValue() + size, j.intValue(),
        //     i.intValue() + size < MAX_SIZE).parent(level);
        // neighbors[2] = fromFaceIJSame(face, i.intValue(), j.intValue() + size,
        //     j.intValue() + size < MAX_SIZE).parent(level);
        // neighbors[3] = fromFaceIJSame(face, i.intValue() - size, j.intValue(),
        //     i.intValue() - size >= 0).parent(level);
        return neighbors;
    };
    /**
     * Return the neighbors of closest vertex to this cell at the given level, by
     * appending them to "output". Normally there are four neighbors, but the
     * closest vertex may only have three neighbors if it is one of the 8 cube
     * vertices.
     *
     * Requires: level < this.evel(), so that we can determine which vertex is
     * closest (in particular, level == MAX_LEVEL is not allowed).
     */
    S2CellId.prototype.getVertexNeighbors = function (level) {
        // "level" must be strictly less than this cell's level so that we can
        // determine which vertex this cell is closest to.
        // assert (level < this.level());
        var i = new MutableInteger(0);
        var j = new MutableInteger(0);
        var face = this.toFaceIJOrientation(i, j, null);
        // Determine the i- and j-offsets to the closest neighboring cell in each
        // direction. This involves looking at the next bit of "i" and "j" to
        // determine which quadrant of this->parent(level) this cell lies in.
        var halfsize = 1 << (S2CellId.MAX_LEVEL - (level + 1));
        var size = halfsize << 1;
        var isame, jsame;
        var ioffset, joffset;
        if ((i.val & halfsize) != 0) {
            ioffset = size;
            isame = (i.val + size) < S2CellId.MAX_SIZE;
        }
        else {
            ioffset = -size;
            isame = (i.val - size) >= 0;
        }
        if ((j.val & halfsize) != 0) {
            joffset = size;
            jsame = (j.val + size) < S2CellId.MAX_SIZE;
        }
        else {
            joffset = -size;
            jsame = (j.val - size) >= 0;
        }
        var toRet = [];
        toRet.push(this.parentL(level));
        toRet.push(S2CellId
            .fromFaceIJSame(face, i.val + ioffset, j.val, isame)
            .parentL(level));
        // output
        //     .add(fromFaceIJSame(face, i.intValue() + ioffset, j.intValue(), isame)
        //         .parent(level));
        toRet.push(S2CellId
            .fromFaceIJSame(face, i.val, j.val + joffset, jsame)
            .parentL(level));
        // output
        //     .add(fromFaceIJSame(face, i.intValue(), j.intValue() + joffset, jsame)
        //         .parent(level));
        // If i- and j- edge neighbors are *both* on a different face, then this
        // vertex only has three neighbors (it is one of the 8 cube vertices).
        if (isame || jsame) {
            toRet.push(S2CellId.fromFaceIJSame(face, i.val + ioffset, j.val + joffset, isame && jsame).parentL(level));
            // output.add(fromFaceIJSame(face, i.intValue() + ioffset,
            //     j.intValue() + joffset, isame && jsame).parent(level));
        }
        return toRet;
    };
    /**
     * Append all neighbors of this cell at the given level to "output". Two cells
     * X and Y are neighbors if their boundaries intersect but their interiors do
     * not. In particular, two cells that intersect at a single point are
     * neighbors.
     *
     * Requires: nbr_level >= this->level(). Note that for cells adjacent to a
     * face vertex, the same neighbor may be appended more than once.
     */
    S2CellId.prototype.getAllNeighbors = function (nbrLevel) {
        var i = new MutableInteger(0);
        var j = new MutableInteger(0);
        var face = this.toFaceIJOrientation(i, j, null);
        // Find the coordinates of the lower left-hand leaf cell. We need to
        // normalize (i,j) to a known position within the cell because nbr_level
        // may be larger than this cell's level.
        var size = 1 << (S2CellId.MAX_LEVEL - this.level());
        i.val = i.val & -size;
        j.val = j.val & -size;
        var nbrSize = 1 << (S2CellId.MAX_LEVEL - nbrLevel);
        // assert (nbrSize <= size);
        var output = [];
        // We compute the N-S, E-W, and diagonal neighbors in one pass.
        // The loop test is at the end of the loop to avoid 32-bit overflow.
        for (var k = -nbrSize;; k += nbrSize) {
            var sameFace = void 0;
            if (k < 0) {
                sameFace = (j.val + k >= 0);
            }
            else if (k >= size) {
                sameFace = (j.val + k < S2CellId.MAX_SIZE);
            }
            else {
                sameFace = true;
                // North and South neighbors.
                output.push(S2CellId.fromFaceIJSame(face, i.val + k, j.val - nbrSize, j.val - size >= 0).parentL(nbrLevel));
                output.push(S2CellId.fromFaceIJSame(face, i.val + k, j.val + size, j.val + size < S2CellId.MAX_SIZE).parentL(nbrLevel));
            }
            // East, West, and Diagonal neighbors.
            output.push(S2CellId.fromFaceIJSame(face, i.val - nbrSize, j.val + k, sameFace && i.val - size >= 0).parentL(nbrLevel));
            output.push(S2CellId.fromFaceIJSame(face, i.val + size, j.val + k, sameFace && i.val + size < S2CellId.MAX_SIZE).parentL(nbrLevel));
            if (k >= size) {
                break;
            }
        }
        return output;
    };
    // ///////////////////////////////////////////////////////////////////
    // Low-level methods.
    /**
     * Return a leaf cell given its cube face (range 0..5) and i- and
     * j-coordinates (see s2.h).
     */
    S2CellId.fromFaceIJ = function (face, i, j) {
        // Optimization notes:
        // - Non-overlapping bit fields can be combined with either "+" or "|".
        // Generally "+" seems to produce better code, but not always.
        // gcc doesn't have very good code generation for 64-bit operations.
        // We optimize this by computing the result as two 32-bit integers
        // and combining them at the end. Declaring the result as an array
        // rather than local variables helps the compiler to do a better job
        // of register allocation as well. Note that the two 32-bits halves
        // get shifted one bit to the left when they are combined.
        var faceL = new Long(face);
        var n = [new Long(0), faceL.shiftLeft(S2CellId.POS_BITS - 33)];
        // Alternating faces have opposite Hilbert curve orientations; this
        // is necessary in order for all faces to have a right-handed
        // coordinate system.
        var bits = faceL.and(S2CellId.SWAP_MASK);
        // Each iteration maps 4 bits of "i" and "j" into 8 bits of the Hilbert
        // curve position. The lookup table transforms a 10-bit key of the form
        // "iiiijjjjoo" to a 10-bit value of the form "ppppppppoo", where the
        // letters [ijpo] denote bits of "i", "j", Hilbert curve position, and
        // Hilbert curve orientation respectively.
        for (var k = 7; k >= 0; --k) {
            bits = S2CellId.getBits(n, i, j, k, bits);
        }
        // S2CellId s = new S2CellId((((n[1] << 32) + n[0]) << 1) + 1);
        return new S2CellId(n[1].shiftLeft(32)
            .add(n[0])
            .shiftLeft(1)
            .add(1));
    };
    S2CellId.getBits = function (n, i, j, k, bits) {
        var mask = new Long(1).shiftLeft(S2CellId.LOOKUP_BITS).sub(1);
        bits = bits.add(new Long(i)
            .shiftRight(k * S2CellId.LOOKUP_BITS)
            .and(mask)
            .shiftLeft(S2CellId.LOOKUP_BITS + 2));
        // bits += (((i >> (k * LOOKUP_BITS)) & mask) << (LOOKUP_BITS + 2));
        bits = bits.add(new Long(j)
            .shiftRight(k * S2CellId.LOOKUP_BITS)
            .and(mask)
            .shiftLeft(2));
        // bits += (((j >> (k * LOOKUP_BITS)) & mask) << 2);
        bits = S2CellId.LOOKUP_POS[bits.toNumber()];
        n[k >> 2] = n[k >> 2].or(bits.shiftRight(2).shiftLeft((k & 3) * 2 * S2CellId.LOOKUP_BITS));
        // n[k >> 2] |= ((((long) bits) >> 2) << ((k & 3) * 2 * LOOKUP_BITS));
        return bits.and(S2CellId.SWAP_MASK | S2CellId.INVERT_MASK);
    };
    /**
     * Return the i- or j-index of the leaf cell containing the given s- or
     * t-value.
     */
    S2CellId.stToIJ = function (_s) {
        // Converting from floating-point to integers via static_cast is very slow
        // on Intel processors because it requires changing the rounding mode.
        // Rounding to the nearest integer using FastIntRound() is much faster.
        var s = S2.toDecimal(_s);
        var m = S2.toDecimal(S2CellId.MAX_SIZE).dividedBy(2); // scaling multiplier
        return decimal.Decimal.max(0, decimal.Decimal.min(m.times(2).minus(1), decimal.Decimal.round(m.times(s).plus(m.minus(0.5))))).toNumber();
        // return Math.max(0,  Math.min(2 * m - 1, Math.round(m * s + (m - 0.5))));
        // return (int) Math.max(0, Math.min(2 * m - 1, Math.round(m * s + (m - 0.5))));
    };
    /**
     * Given (i, j) coordinates that may be out of bounds, normalize them by
     * returning the corresponding neighbor cell on an adjacent face.
     */
    S2CellId.fromFaceIJWrap = function (face, i, j) {
        // Convert i and j to the coordinates of a leaf cell just beyond the
        // boundary of this face. This prevents 32-bit overflow in the case
        // of finding the neighbors of a face cell, and also means that we
        // don't need to worry about the distinction between (s,t) and (u,v).
        i = Math.max(-1, Math.min(S2CellId.MAX_SIZE, i));
        j = Math.max(-1, Math.min(S2CellId.MAX_SIZE, j));
        // Find the (s,t) coordinates corresponding to (i,j). At least one
        // of these coordinates will be just outside the range [0, 1].
        var kScale = S2.toDecimal(1.0).dividedBy(S2CellId.MAX_SIZE);
        var s = kScale.times(new Long(i).shiftLeft(1).add(1).sub(S2CellId.MAX_SIZE).toInt());
        var t = kScale.times(new Long(j).shiftLeft(1).add(1).sub(S2CellId.MAX_SIZE).toInt());
        // Find the leaf cell coordinates on the adjacent face, and convert
        // them to a cell id at the appropriate level.
        var p = new R2Vector(s, t).toPoint(face);
        face = p.toFace();
        // face = S2Projections.xyzToFace(p);
        var st = p.toR2Vector(face);
        // R2Vector st = S2Projections.validFaceXyzToUv(face, p);
        return S2CellId.fromFaceIJ(face, S2CellId.stToIJ(st.x), S2CellId.stToIJ(st.y));
    };
    /**
     * Public helper function that calls FromFaceIJ if sameFace is true, or
     * FromFaceIJWrap if sameFace is false.
     */
    S2CellId.fromFaceIJSame = function (face, i, j, sameFace) {
        if (sameFace) {
            return S2CellId.fromFaceIJ(face, i, j);
        }
        else {
            return S2CellId.fromFaceIJWrap(face, i, j);
        }
    };
    /**
     * Returns true if x1 < x2, when both values are treated as unsigned.
     */
    S2CellId.unsignedLongLessThan = function (x1, x2) {
        return x1.toUnsigned().lessThan(x2.toUnsigned());
        // return (x1 + Long.MIN_VALUE) < (x2 + Long.MIN_VALUE);
    };
    /**
     * Returns true if x1 > x2, when both values are treated as unsigned.
     */
    S2CellId.unsignedLongGreaterThan = function (x1, x2) {
        return x1.toUnsigned().greaterThan(x2.toUnsigned());
        // return (x1 + Long.MIN_VALUE) > (x2 + Long.MIN_VALUE);
    };
    S2CellId.prototype.lessThan = function (x) {
        return S2CellId.unsignedLongLessThan(this.id, x.id);
    };
    S2CellId.prototype.greaterThan = function (x) {
        return S2CellId.unsignedLongGreaterThan(this.id, x.id);
    };
    S2CellId.prototype.lessOrEquals = function (x) {
        return S2CellId.unsignedLongLessThan(this.id, x.id) || this.id.equals(x.id);
    };
    S2CellId.prototype.greaterOrEquals = function (x) {
        return S2CellId.unsignedLongGreaterThan(this.id, x.id) || this.id.equals(x.id);
    };
    S2CellId.prototype.toString = function () {
        return "(face=" + this.face + ", pos=" + this.pos().toString(16) + ", level="
            + this.level() + ")";
    };
    S2CellId.prototype.compareTo = function (that) {
        return S2CellId.unsignedLongLessThan(this.id, that.id) ? -1 :
            S2CellId.unsignedLongGreaterThan(this.id, that.id) ? 1 : 0;
    };
    S2CellId.prototype.equals = function (that) {
        return this.compareTo(that) === 0;
    };
    /**
     * Returns the position of the id within the given list or a negative value with
     * the position of the index wher eit should be entered if the id was present
     */
    S2CellId.binarySearch = function (ids, _id, low) {
        if (low === void 0) { low = 0; }
        var id;
        if (_id instanceof S2CellId) {
            id = _id;
        }
        else if (_id instanceof Long) {
            id = new S2CellId(_id);
        }
        var high = ids.length - 1;
        while (low <= high) {
            var mid = (low + high) >>> 1;
            var midVal = ids[mid];
            var cmp = midVal.compareTo(id);
            if (cmp < 0)
                low = mid + 1;
            else if (cmp > 0)
                high = mid - 1;
            else
                return mid; // key found
        }
        return -(low + 1); // key not found
    };
    S2CellId.indexedBinarySearch = function (ids, id, low) {
        if (low === void 0) { low = 0; }
        var toRet = this.binarySearch(ids, id, low);
        if (toRet >= 0) {
            return toRet;
        }
        else {
            return -(toRet + 1);
        }
    };
    // Although only 60 bits are needed to represent the index of a leaf
    // cell, we need an extra bit in order to represent the position of
    // the center of the leaf cell along the Hilbert curve.
    S2CellId.FACE_BITS = 3;
    S2CellId.NUM_FACES = 6;
    S2CellId.MAX_LEVEL = 30; // Valid levels: 0..MAX_LEVEL
    S2CellId.POS_BITS = 2 * S2CellId.MAX_LEVEL + 1;
    S2CellId.MAX_SIZE = 1 << S2CellId.MAX_LEVEL;
    //
    // calculated as 0xffffffffffffffff / radix
    S2CellId.maxValueDivs = [new Long(0), new Long(0),
        parseHex('9223372036854775807'), parseHex('6148914691236517205'), parseHex('4611686018427387903'),
        parseHex('3689348814741910323'), parseHex('3074457345618258602'), parseHex('2635249153387078802'),
        parseHex('2305843009213693951'), parseHex('2049638230412172401'), parseHex('1844674407370955161'),
        parseHex('1676976733973595601'), parseHex('1537228672809129301'), parseHex('1418980313362273201'),
        parseHex('1317624576693539401'), parseHex('1229782938247303441'), parseHex('1152921504606846975'),
        parseHex('1085102592571150095'), parseHex('1024819115206086200'), parseHex('970881267037344821'),
        parseHex('922337203685477580'), parseHex('878416384462359600'), parseHex('838488366986797800'),
        parseHex('802032351030850070'), parseHex('768614336404564650'), parseHex('737869762948382064'),
        parseHex('709490156681136600'), parseHex('683212743470724133'), parseHex('658812288346769700'),
        parseHex('636094623231363848'), parseHex('614891469123651720'), parseHex('595056260442243600'),
        parseHex('576460752303423487'), parseHex('558992244657865200'), parseHex('542551296285575047'),
        parseHex('527049830677415760'), parseHex('512409557603043100')]; // 35-36
    // calculated as 0xffffffffffffffff % radix
    S2CellId.maxValueMods = [0, 0,
        1, 0, 3, 0, 3, 1, 7, 6, 5, 4, 3, 2, 1, 0, 15, 0, 15, 16, 15, 15,
        15, 5, 15, 15, 15, 24, 15, 23, 15, 15, 31, 15, 17, 15, 15]; // 22-36
    // Constant related to unsigned long's
    // '18446744073709551615'
    // Long.fromString('0xffffffffffffffff', true, 16).toString()
    // new Decimal(2).pow(64).sub(1);
    S2CellId.MAX_UNSIGNED = Long.fromString('0xffffffffffffffff', true, 16);
    // The following lookup tables are used to convert efficiently between an
    // (i,j) cell index and the corresponding position along the Hilbert curve.
    // "lookup_pos" maps 4 bits of "i", 4 bits of "j", and 2 bits representing the
    // orientation of the current cell into 8 bits representing the order in which
    // that subcell is visited by the Hilbert curve, plus 2 bits indicating the
    // new orientation of the Hilbert curve within that subcell. (Cell
    // orientations are represented as combination of kSwapMask and kInvertMask.)
    //
    // "lookup_ij" is an inverted table used for mapping in the opposite
    // direction.
    //
    // We also experimented with looking up 16 bits at a time (14 bits of position
    // plus 2 of orientation) but found that smaller lookup tables gave better
    // performance. (2KB fits easily in the primary cache.)
    // Values for these constants are *declared* in the *.h file. Even though
    // the declaration specifies a value for the constant, that declaration
    // is not a *definition* of storage for the value. Because the values are
    // supplied in the declaration, we don't need the values here. Failing to
    // define storage causes link errors for any code that tries to take the
    // address of one of these values.
    S2CellId.LOOKUP_BITS = 4;
    S2CellId.SWAP_MASK = 0x01;
    S2CellId.INVERT_MASK = 0x02;
    S2CellId.LOOKUP_POS = [];
    S2CellId.LOOKUP_IJ = [];
    /**
     * This is the offset required to wrap around from the beginning of the
     * Hilbert curve to the end or vice versa; see next_wrap() and prev_wrap().
     */
    S2CellId.WRAP_OFFSET = new Long(S2CellId.NUM_FACES).shiftLeft(S2CellId.POS_BITS);
    return S2CellId;
}());
function initLookupCell(level, i, j, origOrientation, pos, orientation) {
    if (level == S2CellId.LOOKUP_BITS) {
        var ij = (i << S2CellId.LOOKUP_BITS) + j;
        S2CellId.LOOKUP_POS[(ij << 2) + origOrientation] = pos.shiftLeft(2).add(orientation);
        S2CellId.LOOKUP_IJ[pos.shiftLeft(2).add(origOrientation).toNumber()] = (ij << 2) + orientation;
        // new Long((ij << 2)).add(orientation);
    }
    else {
        level++;
        i <<= 1;
        j <<= 1;
        pos = pos.shiftLeft(2);
        // Initialize each sub-cell recursively.
        for (var subPos = 0; subPos < 4; subPos++) {
            var ij = S2.POS_TO_IJ[orientation][subPos];
            var orientationMask = S2.POS_TO_ORIENTATION[subPos];
            initLookupCell(level, i + (ij >>> 1), j + (ij & 1), origOrientation, pos.add(subPos), orientation ^ orientationMask);
        }
    }
}
initLookupCell(0, 0, 0, 0, new Long(0), 0);
initLookupCell(0, 0, 0, S2.SWAP_MASK, new Long(0), S2.SWAP_MASK);
initLookupCell(0, 0, 0, S2.INVERT_MASK, new Long(0), S2.INVERT_MASK);
initLookupCell(0, 0, 0, S2.SWAP_MASK | S2.INVERT_MASK, new Long(0), S2.SWAP_MASK | S2.INVERT_MASK);

/*
 * Copyright 2005 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
(function (Projections) {
    Projections[Projections["S2_LINEAR_PROJECTION"] = 0] = "S2_LINEAR_PROJECTION";
    Projections[Projections["S2_TAN_PROJECTION"] = 1] = "S2_TAN_PROJECTION";
    Projections[Projections["S2_QUADRATIC_PROJECTION"] = 2] = "S2_QUADRATIC_PROJECTION";
})(exports.Projections || (exports.Projections = {}));
var S2Projections = /** @class */ (function () {
    function S2Projections() {
    }
    S2Projections.getUNorm = function (face, u) {
        switch (face) {
            case 0:
                return new S2Point(u, -1, 0);
            case 1:
                return new S2Point(1, u, 0);
            case 2:
                return new S2Point(1, 0, u);
            case 3:
                return new S2Point(-u, 0, 1);
            case 4:
                return new S2Point(0, -u, 1);
            default:
                return new S2Point(0, -1, -u);
        }
    };
    S2Projections.getVNorm = function (face, v) {
        switch (face) {
            case 0:
                return new S2Point(-v, 0, 1);
            case 1:
                return new S2Point(0, -v, 1);
            case 2:
                return new S2Point(0, -1, -v);
            case 3:
                return new S2Point(v, -1, 0);
            case 4:
                return new S2Point(1, v, 0);
            default:
                return new S2Point(1, 0, v);
        }
    };
    S2Projections.getUAxis = function (face) {
        switch (face) {
            case 0:
                return new S2Point(0, 1, 0);
            case 1:
                return new S2Point(-1, 0, 0);
            case 2:
                return new S2Point(-1, 0, 0);
            case 3:
                return new S2Point(0, 0, -1);
            case 4:
                return new S2Point(0, 0, -1);
            default:
                return new S2Point(0, 1, 0);
        }
    };
    S2Projections.getVAxis = function (face) {
        switch (face) {
            case 0:
                return new S2Point(0, 0, 1);
            case 1:
                return new S2Point(0, 0, 1);
            case 2:
                return new S2Point(0, -1, 0);
            case 3:
                return new S2Point(0, -1, 0);
            case 4:
                return new S2Point(1, 0, 0);
            default:
                return new S2Point(1, 0, 0);
        }
    };
    S2Projections.faceUvToXyz = function (face, u, v) {
        return new R2Vector(u, v).toPoint(face);
    };
    S2Projections.MIN_WIDTH = new S2Metric$$1(1, S2.M_SQRT2 / 3);
    S2Projections.AVG_AREA = new S2Metric$$1(2, S2.M_PI / 6); // 0.524)
    return S2Projections;
}());

var S2Cell = /** @class */ (function () {
    function S2Cell(cellID) {
        this.cellID = cellID;
        this._uv = [];
        this._uv.push([]);
        this._uv.push([]);
        this.init(cellID);
    }
    Object.defineProperty(S2Cell.prototype, "id", {
        get: function () {
            return this.cellID;
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(S2Cell.prototype, "face", {
        get: function () {
            return this._face;
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(S2Cell.prototype, "level", {
        get: function () {
            return this._level;
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(S2Cell.prototype, "orientation", {
        get: function () {
            return this._orientation;
        },
        enumerable: true,
        configurable: true
    });
    // This is a static method in order to provide named parameters.
    S2Cell.fromFacePosLevel = function (face, pos, level) {
        return new S2Cell(S2CellId.fromFacePosLevel(face, new Long(pos), level));
    };
    // Convenience methods.
    S2Cell.fromPoint = function (p) {
        return new S2Cell(S2CellId.fromPoint(p));
    };
    S2Cell.fromLatLng = function (ll) {
        return new S2Cell(S2CellId.fromPoint(ll.toPoint()));
    };
    S2Cell.prototype.isLeaf = function () {
        return this.level == S2CellId.MAX_LEVEL;
    };
    S2Cell.prototype.getVertex = function (k) {
        return S2Point.normalize(this.getVertexRaw(k));
    };
    /**
     * Return the k-th vertex of the cell (k = 0,1,2,3). Vertices are returned in
     * CCW order. The points returned by GetVertexRaw are not necessarily unit
     * length.
     */
    S2Cell.prototype.getVertexRaw = function (k) {
        // Vertices are returned in the order SW, SE, NE, NW.
        return new R2Vector(this._uv[0][(k >> 1) ^ (k & 1)], this._uv[1][k >> 1])
            .toPoint(this.face);
        // return S2Projections.faceUvToXyz(this.face, );
    };
    S2Cell.prototype.getEdge = function (k) {
        return S2Point.normalize(this.getEdgeRaw(k));
    };
    S2Cell.prototype.getEdgeRaw = function (k) {
        switch (k) {
            case 0:
                return S2Projections.getVNorm(this.face, this._uv[1][0]); // South
            case 1:
                return S2Projections.getUNorm(this.face, this._uv[0][1]); // East
            case 2:
                return S2Point.neg(S2Projections.getVNorm(this.face, this._uv[1][1])); // North
            default:
                return S2Point.neg(S2Projections.getUNorm(this.face, this._uv[0][0])); // West
        }
    };
    /**
     * Return the inward-facing normal of the great circle passing through the
     * edge from vertex k to vertex k+1 (mod 4). The normals returned by
     * GetEdgeRaw are not necessarily unit length.
     *
     *  If this is not a leaf cell, set children[0..3] to the four children of
     * this cell (in traversal order) and return true. Otherwise returns false.
     * This method is equivalent to the following:
     *
     *  for (pos=0, id=child_begin(); id != child_end(); id = id.next(), ++pos)
     * children[i] = S2Cell(id);
     *
     * except that it is more than two times faster.
     */
    S2Cell.prototype.subdivide = function () {
        // This function is equivalent to just iterating over the child cell ids
        // and calling the S2Cell constructor, but it is about 2.5 times faster.
        if (this.isLeaf()) {
            return null;
        }
        // Compute the cell midpoint in uv-space.
        // const uvMid = this.getCenterUV();
        var children = new Array(4);
        // Create four children with the appropriate bounds.
        var id = this.cellID.childBegin();
        for (var pos = 0; pos < 4; ++pos, id = id.next()) {
            children[pos] = new S2Cell(id);
            // S2Cell child = children[pos];
            // child.face = this.face;
            // child.level = (byte) (this.level + 1);
            // child.orientation = (byte) (this.orientation ^ S2.posToOrientation(pos));
            // child.cellId = id;
            // int ij = S2.posToIJ(this.orientation, pos);
            // for (let d = 0; d < 2; ++d) {
            //   // The dimension 0 index (i/u) is in bit 1 of ij.
            //   int m = 1 - ((ij >> (1 - d)) & 1);
            //   child._uv[d][m] = uvMid.get(d);
            //   child._uv[d][1 - m] = this._uv[d][1 - m];
            // }
        }
        return children;
    };
    /**
     * Return the direction vector corresponding to the center in (s,t)-space of
     * the given cell. This is the point at which the cell is divided into four
     * subcells; it is not necessarily the centroid of the cell in (u,v)-space or
     * (x,y,z)-space. The point returned by GetCenterRaw is not necessarily unit
     * length.
     */
    S2Cell.prototype.getCenter = function () {
        return S2Point.normalize(this.getCenterRaw());
    };
    S2Cell.prototype.getCenterRaw = function () {
        return this.cellID.toPointRaw();
    };
    /**
     * Return the center of the cell in (u,v) coordinates (see {@code
     * S2Projections}). Note that the center of the cell is defined as the point
     * at which it is recursively subdivided into four children; in general, it is
     * not at the midpoint of the (u,v) rectangle covered by the cell
     */
    S2Cell.prototype.getCenterUV = function () {
        var i = new MutableInteger(0);
        var j = new MutableInteger(0);
        this.cellID.toFaceIJOrientation(i, j, null);
        var cellSize = 1 << (S2CellId.MAX_LEVEL - this.level);
        // TODO(dbeaumont): Figure out a better naming of the variables here (and elsewhere).
        var si = (i.val & -cellSize) * 2 + cellSize - S2Cell.MAX_CELL_SIZE;
        var x = R2Vector.singleStTOUV(S2.toDecimal(1).dividedBy(S2Cell.MAX_CELL_SIZE).times(si));
        // let x = S2Projections.stToUV((1.0 / S2Cell.MAX_CELL_SIZE) * si);
        var sj = (j.val & -cellSize) * 2 + cellSize - S2Cell.MAX_CELL_SIZE;
        var y = R2Vector.singleStTOUV(S2.toDecimal(1).dividedBy(S2Cell.MAX_CELL_SIZE).times(sj));
        // double y = S2Projections.stToUV((1.0 / S2Cell.MAX_CELL_SIZE) * sj);
        return new R2Vector(x, y);
    };
    /**
     * Return the average area of cells at this level. This is accurate to within
     * a factor of 1.7 (for S2_QUADRATIC_PROJECTION) and is extremely cheap to
     * compute.
     */
    S2Cell.averageArea = function (level) {
        return S2Projections.AVG_AREA.getValue(level);
    };
    /**
     * Return the average area of cells at this level. This is accurate to within
     * a factor of 1.7 (for S2_QUADRATIC_PROJECTION) and is extremely cheap to
     * compute.
     */
    S2Cell.prototype.averageArea = function () {
        return S2Projections.AVG_AREA.getValue(this.level);
    };
    /**
     * Return the approximate area of this cell. This method is accurate to within
     * 3% percent for all cell sizes and accurate to within 0.1% for cells at
     * level 5 or higher (i.e. 300km square or smaller). It is moderately cheap to
     * compute.
     */
    S2Cell.prototype.approxArea = function () {
        // All cells at the first two levels have the same area.
        if (this.level < 2) {
            return this.averageArea();
        }
        // First, compute the approximate area of the cell when projected
        // perpendicular to its normal. The cross product of its diagonals gives
        // the normal, and the length of the normal is twice the projected area.
        var flatArea = S2Point.crossProd(S2Point.sub(this.getVertex(2), this.getVertex(0)), S2Point.sub(this.getVertex(3), this.getVertex(1))).norm().times(0.5);
        // double flatArea = 0.5 * S2Point.crossProd(
        //         S2Point.sub(getVertex(2), getVertex(0)), S2Point.sub(getVertex(3), getVertex(1))).norm();
        // Now, compensate for the curvature of the cell surface by pretending
        // that the cell is shaped like a spherical cap. The ratio of the
        // area of a spherical cap to the area of its projected disc turns out
        // to be 2 / (1 + sqrt(1 - r*r)) where "r" is the radius of the disc.
        // For example, when r=0 the ratio is 1, and when r=1 the ratio is 2.
        // Here we set Pi*r*r == flat_area to find the equivalent disc.
        return flatArea
            .times(2)
            .dividedBy(decimal.Decimal.min(flatArea.times(S2.M_1_PI), 1)
            .neg()
            .plus(1)
            .sqrt()
            .plus(1)).toNumber();
    };
    //
    // /**
    //  * Return the area of this cell as accurately as possible. This method is more
    //  * expensive but it is accurate to 6 digits of precision even for leaf cells
    //  * (whose area is approximately 1e-18).
    //  */
    S2Cell.prototype.exactArea = function () {
        var v0 = this.getVertex(0);
        var v1 = this.getVertex(1);
        var v2 = this.getVertex(2);
        var v3 = this.getVertex(3);
        return S2.area(v0, v1, v2).plus(S2.area(v0, v2, v3));
    };
    // //////////////////////////////////////////////////////////////////////
    // S2Region interface (see {@code S2Region} for details):
    S2Cell.prototype.getCapBound = function () {
        // Use the cell center in (u,v)-space as the cap axis. This vector is
        // very close to GetCenter() and faster to compute. Neither one of these
        // vectors yields the bounding cap with minimal surface area, but they
        // are both pretty close.
        //
        // It's possible to show that the two vertices that are furthest from
        // the (u,v)-origin never determine the maximum cap size (this is a
        // possible future optimization).
        var u = this._uv[0][0].plus(this._uv[0][1]).times(0.5);
        var v = this._uv[1][0].plus(this._uv[1][1]).times(0.5);
        var cap = new S2Cap(S2Point.normalize(S2Projections.faceUvToXyz(this.face, u, v)), 0);
        for (var k = 0; k < 4; ++k) {
            cap = cap.addPoint(this.getVertex(k));
        }
        return cap;
    };
    // 35.26 degrees
    S2Cell.prototype.getRectBound = function () {
        if (this.level > 0) {
            // Except for cells at level 0, the latitude and longitude extremes are
            // attained at the vertices. Furthermore, the latitude range is
            // determined by one pair of diagonally opposite vertices and the
            // longitude range is determined by the other pair.
            //
            // We first determine which corner (i,j) of the cell has the largest
            // absolute latitude. To maximize latitude, we want to find the point in
            // the cell that has the largest absolute z-coordinate and the smallest
            // absolute x- and y-coordinates. To do this we look at each coordinate
            // (u and v), and determine whether we want to minimize or maximize that
            // coordinate based on the axis direction and the cell's (u,v) quadrant.
            var u = this._uv[0][0].plus(this._uv[0][1]);
            var v = this._uv[1][0].plus(this._uv[1][1]);
            var i = S2Projections.getUAxis(this.face).z.eq(0) ? (u.lt(0) ? 1 : 0) : (u.gt(0) ? 1 : 0);
            var j = S2Projections.getVAxis(this.face).z.eq(0) ? (v.lt(0) ? 1 : 0) : (v.gt(0) ? 1 : 0);
            var lat = R1Interval.fromPointPair(this.getLatitude(i, j), this.getLatitude(1 - i, 1 - j));
            lat = lat.expanded(S2Cell.MAX_ERROR).intersection(S2LatLngRect.fullLat());
            if (lat.lo.eq(-S2.M_PI_2) || lat.hi.eq(S2.M_PI_2)) {
                return new S2LatLngRect(lat, S1Interval.full());
            }
            var lng = S1Interval.fromPointPair(this.getLongitude(i, 1 - j), this.getLongitude(1 - i, j));
            return new S2LatLngRect(lat, lng.expanded(S2Cell.MAX_ERROR));
        }
        // The face centers are the +X, +Y, +Z, -X, -Y, -Z axes in that order.
        // assert (S2Projections.getNorm(face).get(face % 3) == ((face < 3) ? 1 : -1));
        switch (this.face) {
            case 0:
                return new S2LatLngRect(new R1Interval(-S2.M_PI_4, S2.M_PI_4), new S1Interval(-S2.M_PI_4, S2.M_PI_4));
            case 1:
                return new S2LatLngRect(new R1Interval(-S2.M_PI_4, S2.M_PI_4), new S1Interval(S2.M_PI_4, 3 * S2.M_PI_4));
            case 2:
                return new S2LatLngRect(new R1Interval(S2Cell.POLE_MIN_LAT, S2.M_PI_2), new S1Interval(-S2.M_PI, S2.M_PI));
            case 3:
                return new S2LatLngRect(new R1Interval(-S2.M_PI_4, S2.M_PI_4), new S1Interval(3 * S2.M_PI_4, -3 * S2.M_PI_4));
            case 4:
                return new S2LatLngRect(new R1Interval(-S2.M_PI_4, S2.M_PI_4), new S1Interval(-3 * S2.M_PI_4, -S2.M_PI_4));
            default:
                return new S2LatLngRect(new R1Interval(-S2.M_PI_2, -S2Cell.POLE_MIN_LAT), new S1Interval(-S2.M_PI, S2.M_PI));
        }
    };
    S2Cell.prototype.mayIntersect = function (cell) {
        return this.cellID.intersects(cell.cellID);
    };
    S2Cell.prototype.contains = function (p) {
        // We can't just call XYZtoFaceUV, because for points that lie on the
        // boundary between two faces (i.e. u or v is +1/-1) we need to return
        // true for both adjacent cells.
        var uvPoint = p.toR2Vector(this.face);
        // S2Projections.faceXyzToUv(this.face, p);
        if (uvPoint == null) {
            return false;
        }
        return (uvPoint.x.gte(this._uv[0][0]) && uvPoint.x.lte(this._uv[0][1])
            && uvPoint.y.gte(this._uv[1][0]) && uvPoint.y.lte(this._uv[1][1]));
    };
    // The point 'p' does not need to be normalized.
    S2Cell.prototype.containsC = function (cell) {
        return this.cellID.contains(cell.cellID);
    };
    S2Cell.prototype.init = function (id) {
        this.cellID = id;
        var ij = [];
        var mOrientation = new MutableInteger(0);
        for (var d = 0; d < 2; ++d) {
            ij[d] = new MutableInteger(0);
        }
        this._face = id.toFaceIJOrientation(ij[0], ij[1], mOrientation);
        this._orientation = mOrientation.val; // Compress int to a byte.
        this._level = id.level();
        var cellSize = 1 << (S2CellId.MAX_LEVEL - this.level);
        for (var d = 0; d < 2; ++d) {
            // Compute the cell bounds in scaled (i,j) coordinates.
            var sijLo = (ij[d].val & -cellSize) * 2 - S2Cell.MAX_CELL_SIZE;
            var sijHi = sijLo + cellSize * 2;
            var s = S2.toDecimal(1).dividedBy(S2Cell.MAX_CELL_SIZE);
            this._uv[d][0] = R2Vector.singleStTOUV(s.times(sijLo));
            //S2Projections.stToUV((1.0 / S2Cell.MAX_CELL_SIZE) * sijLo);
            this._uv[d][1] = R2Vector.singleStTOUV(s.times(sijHi));
            //S2Projections.stToUV((1.0 / S2Cell.MAX_CELL_SIZE) * sijHi);
        }
    };
    // Internal method that does the actual work in the constructors.
    S2Cell.prototype.getLatitude = function (i, j) {
        var p = S2Projections.faceUvToXyz(this.face, this._uv[0][i], this._uv[1][j]);
        return decimal.Decimal.atan2(p.z, p.x.pow(2).plus(p.y.pow(2))
            .sqrt());
        // return Math.atan2(p.z, Math.sqrt(p.x * p.x + p.y * p.y));
    };
    S2Cell.prototype.getLongitude = function (i, j) {
        var p = S2Projections.faceUvToXyz(this.face, this._uv[0][i], this._uv[1][j]);
        return decimal.Decimal.atan2(p.y, p.x);
        // Math.atan2(p.y, p.x);
    };
    // Return the latitude or longitude of the cell vertex given by (i,j),
    // where "i" and "j" are either 0 or 1.
    S2Cell.prototype.toString = function () {
        return "[" + this._face + ", " + this._level + ", " + this._orientation + ", " + this.cellID.toToken() + "]";
    };
    S2Cell.prototype.toGEOJSON = function () {
        var coords = [this.getVertex(0), this.getVertex(1), this.getVertex(2), this.getVertex(3), this.getVertex(0)]
            .map(function (v) { return S2LatLng.fromPoint(v); })
            .map(function (v) { return ([v.lngDegrees.toNumber(), v.latDegrees.toNumber()]); });
        // const rectJSON = this.getRectBound().toGEOJSON();
        return {
            type: 'Feature',
            geometry: {
                type: 'Polygon',
                coordinates: [coords]
            },
            properties: {},
            title: "Cell: " + this.id.toToken() + " lvl: " + this.level
        };
        // rectJSON.title = `Cell: ${this.id.toToken()}`;
        // return rectJSON;
    };
    S2Cell.MAX_CELL_SIZE = 1 << S2CellId.MAX_LEVEL;
    // We grow the bounds slightly to make sure that the bounding rectangle
    // also contains the normalized versions of the vertices. Note that the
    // maximum result magnitude is Pi, with a floating-point exponent of 1.
    // Therefore adding or subtracting 2**-51 will always change the result.
    S2Cell.MAX_ERROR = S2.toDecimal(1.0).dividedBy(S2.toDecimal(new Long(1).shiftLeft(51).toString()));
    // The 4 cells around the equator extend to +/-45 degrees latitude at the
    // midpoints of their top and bottom edges. The two cells covering the
    // poles extend down to +/-35.26 degrees at their vertices.
    // adding kMaxError (as opposed to the C version) because of asin and atan2
    // roundoff errors
    S2Cell.POLE_MIN_LAT = decimal.Decimal.asin(S2.toDecimal(1.0).dividedBy(3).sqrt()).minus(S2Cell.MAX_ERROR);
    return S2Cell;
}());

/*
 * Copyright 2005 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/**
 * An S2CellUnion is a region consisting of cells of various sizes. Typically a
 * cell union is used to approximate some other shape. There is a tradeoff
 * between the accuracy of the approximation and how many cells are used. Unlike
 * polygons, cells have a fixed hierarchical structure. This makes them more
 * suitable for optimizations based on preprocessing.
 *
 */
var S2CellUnion = /** @class */ (function () {
    function S2CellUnion() {
        /** The CellIds that form the Union */
        this.cellIds = [];
    }
    S2CellUnion.prototype.S2CellUnion = function () {
    };
    /**
     * Populates a cell union with the given S2CellIds or 64-bit cells ids, and
     * then calls Normalize(). The InitSwap() version takes ownership of the
     * vector data without copying and clears the given vector. These methods may
     * be called multiple times.
     */
    S2CellUnion.prototype.initFromIds = function (cellIds) {
        this.initRawIds(cellIds);
        this.normalize();
    };
    S2CellUnion.prototype.initSwap = function (cellIds) {
        this.initRawSwap(cellIds);
        this.normalize();
    };
    S2CellUnion.prototype.initRawCellIds = function (cellIds) {
        this.cellIds = cellIds;
    };
    S2CellUnion.prototype.initRawIds = function (cellIds) {
        var size = cellIds.length;
        this.cellIds = [];
        for (var i = 0; i < size; i++) {
            this.cellIds.push(new S2CellId(cellIds[i]));
        }
    };
    /**
     * Like Init(), but does not call Normalize(). The cell union *must* be
     * normalized before doing any calculations with it, so it is the caller's
     * responsibility to make sure that the input is normalized. This method is
     * useful when converting cell unions to another representation and back.
     * These methods may be called multiple times.
     */
    S2CellUnion.prototype.initRawSwap = function (cellIds) {
        this.cellIds = [].concat(cellIds);
    };
    S2CellUnion.prototype.size = function () {
        return this.cellIds.length;
    };
    /** Convenience methods for accessing the individual cell ids. */
    S2CellUnion.prototype.cellId = function (i) {
        return this.cellIds[i];
    };
    S2CellUnion.prototype.getCellIds = function () {
        return this.cellIds;
    };
    /**
     * Replaces "output" with an expanded version of the cell union where any
     * cells whose level is less than "min_level" or where (level - min_level) is
     * not a multiple of "level_mod" are replaced by their children, until either
     * both of these conditions are satisfied or the maximum level is reached.
     *
     *  This method allows a covering generated by S2RegionCoverer using
     * min_level() or level_mod() constraints to be stored as a normalized cell
     * union (which allows various geometric computations to be done) and then
     * converted back to the original list of cell ids that satisfies the desired
     * constraints.
     */
    S2CellUnion.prototype.denormalize = function (minLevel, levelMod) {
        // assert (minLevel >= 0 && minLevel <= S2CellId.MAX_LEVEL);
        // assert (levelMod >= 1 && levelMod <= 3);
        var output = [];
        for (var i = 0; i < this.cellIds.length; i++) {
            var id = this.cellIds[i];
            var level = id.level();
            var newLevel = Math.max(minLevel, level);
            if (levelMod > 1) {
                // Round up so that (new_level - min_level) is a multiple of level_mod.
                // (Note that S2CellId::kMaxLevel is a multiple of 1, 2, and 3.)
                newLevel += (S2CellId.MAX_LEVEL - (newLevel - minLevel)) % levelMod;
                newLevel = Math.min(S2CellId.MAX_LEVEL, newLevel);
            }
            if (newLevel == level) {
                output.push(id);
            }
            else {
                var end = id.childEndL(newLevel);
                for (var iid = id.childBeginL(newLevel); !iid.equals(end); iid = iid.next()) {
                    output.push(iid);
                }
            }
        }
        return output;
    };
    /**
     * If there are more than "excess" elements of the cell_ids() vector that are
     * allocated but unused, reallocate the array to eliminate the excess space.
     * This reduces memory usage when many cell unions need to be held in memory
     * at once.
     */
    S2CellUnion.prototype.pack = function () {
        throw new Error('useless');
        // this.cellIds.trimToSize();
    };
    S2CellUnion.prototype.containsC = function (cell) {
        return this.containsCell(cell);
    };
    S2CellUnion.prototype.mayIntersectC = function (cell) {
        return this.mayIntersectCell(cell);
    };
    /**
     * Return true if the cell union contains the given cell id. Containment is
     * defined with respect to regions, e.g. a cell contains its 4 children. This
     * is a fast operation (logarithmic in the size of the cell union).
     */
    S2CellUnion.prototype.contains = function (id) {
        // This function requires that Normalize has been called first.
        //
        // This is an exact test. Each cell occupies a linear span of the S2
        // space-filling curve, and the cell id is simply the position at the center
        // of this span. The cell union ids are sorted in increasing order along
        // the space-filling curve. So we simply find the pair of cell ids that
        // surround the given cell id (using binary search). There is containment
        // if and only if one of these two cell ids contains this cell.
        var pos = S2CellId.binarySearch(this.cellIds, id.id);
        if (pos < 0) {
            pos = -pos - 1;
        }
        if (pos < this.cellIds.length && this.cellIds[pos].rangeMin().lessOrEquals(id)) {
            return true;
        }
        return pos != 0 && this.cellIds[pos - 1].rangeMax().greaterOrEquals(id);
    };
    /**
     * Return true if the cell union intersects the given cell id. This is a fast
     * operation (logarithmic in the size of the cell union).
     */
    S2CellUnion.prototype.intersects = function (id) {
        // This function requires that Normalize has been called first.
        // This is an exact test; see the comments for Contains() above.
        var pos = S2CellId.binarySearch(this.cellIds, id.id);
        if (pos < 0) {
            pos = -pos - 1;
        }
        if (pos < this.cellIds.length && this.cellIds[pos].rangeMin().lessOrEquals(id.rangeMax())) {
            return true;
        }
        return pos != 0 && this.cellIds[pos - 1].rangeMax().greaterOrEquals(id.rangeMin());
    };
    S2CellUnion.prototype.containsUnion = function (that) {
        // A divide-and-conquer or alternating-skip-search approach
        // may be significantly faster in both the average and worst case.
        for (var i = 0; i < that.cellIds.length; i++) {
            if (!this.contains(that.cellIds[i])) {
                return false;
            }
        }
        return true;
    };
    /** This is a fast operation (logarithmic in the size of the cell union). */
    S2CellUnion.prototype.containsCell = function (cell) {
        return this.contains(cell.id);
    };
    /**
     * Return true if this cell union contain/intersects the given other cell
     * union.
     */
    S2CellUnion.prototype.intersectsUnion = function (that) {
        // A divide-and-conquer or alternating-skip-search approach
        // may be significantly faster in both the average and worst case.
        for (var i = 0; i < that.cellIds.length; i++) {
            if (!this.intersects(that.cellIds[i])) {
                return false;
            }
        }
        return true;
    };
    S2CellUnion.prototype.getUnion = function (x, y) {
        // assert (x != this && y != this);
        this.cellIds = [].concat(x.cellIds).concat(y.cellIds);
        this.normalize();
    };
    /**
     * Specialized version of GetIntersection() that gets the intersection of a
     * cell union with the given cell id. This can be useful for "splitting" a
     * cell union into chunks.
     */
    S2CellUnion.prototype.getIntersection = function (x, id) {
        // assert (x != this);
        this.cellIds = [];
        if (x.contains(id)) {
            this.cellIds.push(id);
        }
        else {
            var pos = S2CellId.binarySearch(x.cellIds, id.rangeMin().id);
            if (pos < 0) {
                pos = -pos - 1;
            }
            var idmax = id.rangeMax();
            var size = x.cellIds.length;
            while (pos < size && x.cellIds[pos].lessOrEquals(idmax)) {
                this.cellIds.push(x.cellIds[(pos++)]);
            }
        }
    };
    /**
     * Initialize this cell union to the union or intersection of the two given
     * cell unions. Requires: x != this and y != this.
     */
    S2CellUnion.prototype.getIntersectionUU = function (x, y) {
        // assert (x != this && y != this);
        // This is a fairly efficient calculation that uses binary search to skip
        // over sections of both input vectors. It takes constant time if all the
        // cells of "x" come before or after all the cells of "y" in S2CellId order.
        this.cellIds = [];
        var i = 0;
        var j = 0;
        while (i < x.cellIds.length && j < y.cellIds.length) {
            var imin = x.cellId(i).rangeMin();
            var jmin = y.cellId(j).rangeMin();
            if (imin.greaterThan(jmin)) {
                // Either j->contains(*i) or the two cells are disjoint.
                if (x.cellId(i).lessOrEquals(y.cellId(j).rangeMax())) {
                    this.cellIds.push(x.cellId(i++));
                }
                else {
                    // Advance "j" to the first cell possibly contained by *i.
                    j = S2CellId.indexedBinarySearch(y.cellIds, imin, j + 1);
                    // The previous cell *(j-1) may now contain *i.
                    if (x.cellId(i).lessOrEquals(y.cellId(j - 1).rangeMax())) {
                        --j;
                    }
                }
            }
            else if (jmin.greaterThan(imin)) {
                // Identical to the code above with "i" and "j" reversed.
                if (y.cellId(j).lessOrEquals(x.cellId(i).rangeMax())) {
                    this.cellIds.push(y.cellId(j++));
                }
                else {
                    i = S2CellId.indexedBinarySearch(x.cellIds, jmin, i + 1);
                    if (y.cellId(j).lessOrEquals(x.cellId(i - 1).rangeMax())) {
                        --i;
                    }
                }
            }
            else {
                // "i" and "j" have the same range_min(), so one contains the other.
                if (x.cellId(i).lessThan(y.cellId(j))) {
                    this.cellIds.push(x.cellId(i++));
                }
                else {
                    this.cellIds.push(y.cellId(j++));
                }
            }
        }
        // The output is generated in sorted order, and there should not be any
        // cells that can be merged (provided that both inputs were normalized).
        // assert (!normalize());
    };
    /**
     * Expands the cell union such that it contains all cells of the given level
     * that are adjacent to any cell of the original union. Two cells are defined
     * as adjacent if their boundaries have any points in common, i.e. most cells
     * have 8 adjacent cells (not counting the cell itself).
     *
     *  Note that the size of the output is exponential in "level". For example,
     * if level == 20 and the input has a cell at level 10, there will be on the
     * order of 4000 adjacent cells in the output. For most applications the
     * Expand(min_fraction, min_distance) method below is easier to use.
     */
    S2CellUnion.prototype.expand = function (level) {
        var output = [];
        var levelLsb = S2CellId.lowestOnBitForLevel(level);
        var i = this.size() - 1;
        do {
            var id = this.cellId(i);
            if (id.lowestOnBit().lessThan(levelLsb)) {
                id = id.parentL(level);
                // Optimization: skip over any cells contained by this one. This is
                // especially important when very small regions are being expanded.
                while (i > 0 && id.contains(this.cellId(i - 1))) {
                    --i;
                }
            }
            output.push(id);
            output = output.concat(id.getAllNeighbors(level));
        } while (--i >= 0);
        this.initSwap(output);
    };
    /**
     * Expand the cell union such that it contains all points whose distance to
     * the cell union is at most minRadius, but do not use cells that are more
     * than maxLevelDiff levels higher than the largest cell in the input. The
     * second parameter controls the tradeoff between accuracy and output size
     * when a large region is being expanded by a small amount (e.g. expanding
     * Canada by 1km).
     *
     *  For example, if maxLevelDiff == 4, the region will always be expanded by
     * approximately 1/16 the width of its largest cell. Note that in the worst
     * case, the number of cells in the output can be up to 4 * (1 + 2 **
     * maxLevelDiff) times larger than the number of cells in the input.
     */
    S2CellUnion.prototype.expandA = function (minRadius, maxLevelDiff) {
        var minLevel = S2CellId.MAX_LEVEL;
        for (var i = 0; i < this.cellIds.length; i++) {
            var id = this.cellId(i);
            minLevel = Math.min(minLevel, id.level());
        }
        // Find the maximum level such that all cells are at least "min_radius"
        // wide.
        var radiusLevel = S2Projections.MIN_WIDTH.getMaxLevel(minRadius.radians);
        if (radiusLevel == 0 && minRadius.radians.gt(S2Projections.MIN_WIDTH.getValue(0))) {
            // The requested expansion is greater than the width of a face cell.
            // The easiest way to handle this is to expand twice.
            this.expand(0);
        }
        this.expand(Math.min(minLevel + maxLevelDiff, radiusLevel));
    };
    S2CellUnion.prototype.getCapBound = function () {
        // Compute the approximate centroid of the region. This won't produce the
        // bounding cap of minimal area, but it should be close enough.
        if (this.cellIds.length == 0) {
            return S2Cap.empty();
        }
        var centroid = new S2Point(0, 0, 0);
        this.cellIds.forEach(function (id) {
            var area = S2Cell.averageArea(id.level());
            centroid = S2Point.add(centroid, S2Point.mul(id.toPoint(), area));
        });
        if (centroid.equals(new S2Point(0, 0, 0))) {
            centroid = new S2Point(1, 0, 0);
        }
        else {
            centroid = S2Point.normalize(centroid);
        }
        // Use the centroid as the cap axis, and expand the cap angle so that it
        // contains the bounding caps of all the individual cells. Note that it is
        // *not* sufficient to just bound all the cell vertices because the bounding
        // cap may be concave (i.e. cover more than one hemisphere).
        var cap = new S2Cap(centroid, 0);
        this.cellIds.forEach(function (id) {
            cap = cap.addCap(new S2Cell(id).getCapBound());
        });
        return cap;
    };
    S2CellUnion.prototype.getRectBound = function () {
        var bound = S2LatLngRect.empty();
        this.cellIds.forEach(function (id) {
            bound = bound.union(new S2Cell(id).getRectBound());
        });
        return bound;
    };
    /** This is a fast operation (logarithmic in the size of the cell union). */
    S2CellUnion.prototype.mayIntersectCell = function (cell) {
        return this.intersects(cell.id);
    };
    /**
     * The point 'p' does not need to be normalized. This is a fast operation
     * (logarithmic in the size of the cell union).
     */
    S2CellUnion.prototype.containsPoint = function (p) {
        return this.contains(S2CellId.fromPoint(p));
    };
    /**
     * The number of leaf cells covered by the union.
     * This will be no more than 6*2^60 for the whole sphere.
     *
     * @return the number of leaf cells covered by the union
     */
    S2CellUnion.prototype.leafCellsCovered = function () {
        var numLeaves = new Long(0);
        this.cellIds.forEach(function (id) {
            var invertedLevel = S2CellId.MAX_LEVEL - id.level();
            numLeaves = numLeaves
                .add(new Long(1).shiftLeft(invertedLevel << 1));
        });
        return numLeaves;
    };
    /**
     * Approximate this cell union's area by summing the average area of
     * each contained cell's average area, using {@link S2Cell#averageArea()}.
     * This is equivalent to the number of leaves covered, multiplied by
     * the average area of a leaf.
     * Note that {@link S2Cell#averageArea()} does not take into account
     * distortion of cell, and thus may be off by up to a factor of 1.7.
     * NOTE: Since this is proportional to LeafCellsCovered(), it is
     * always better to use the other function if all you care about is
     * the relative average area between objects.
     *
     * @return the sum of the average area of each contained cell's average area
     */
    S2CellUnion.prototype.averageBasedArea = function () {
        return S2.toDecimal(this.leafCellsCovered().toString()).times(S2Projections.AVG_AREA.getValue(S2CellId.MAX_LEVEL)).toNumber();
    };
    /**
     * Calculates this cell union's area by summing the approximate area for each
     * contained cell, using {@link S2Cell#approxArea()}.
     *
     * @return approximate area of the cell union
     */
    S2CellUnion.prototype.approxArea = function () {
        var area = S2.toDecimal(0);
        this.cellIds.forEach(function (id) {
            area = area.plus(new S2Cell(id).approxArea());
        });
        return area.toNumber();
    };
    /**
     * Calculates this cell union's area by summing the exact area for each
     * contained cell, using the {@link S2Cell#exactArea()}.
     *
     * @return the exact area of the cell union
     */
    S2CellUnion.prototype.exactArea = function () {
        var area = S2.toDecimal(0);
        this.cellIds.forEach(function (id) {
            area = area.plus(new S2Cell(id).exactArea());
        });
        return area.toNumber();
    };
    /**
     * Normalizes the cell union by discarding cells that are contained by other
     * cells, replacing groups of 4 child cells by their parent cell whenever
     * possible, and sorting all the cell ids in increasing order. Returns true if
     * the number of cells was reduced.
     *
     *  This method *must* be called before doing any calculations on the cell
     * union, such as Intersects() or Contains().
     *
     * @return true if the normalize operation had any effect on the cell union,
     *         false if the union was already normalized
     */
    S2CellUnion.prototype.normalize = function () {
        // Optimize the representation by looking for cases where all subcells
        // of a parent cell are present.
        var output = [];
        // ArrayList<S2CellId> output = new ArrayList<>(this.cellIds.size());
        // output.ensureCapacity(this.cellIds.size());
        this.cellIds.sort(function (a, b) { return a.compareTo(b); });
        // Collections.sort(this.cellIds);
        this.cellIds.forEach(function (id) {
            var size = output.length;
            // Check whether this cell is contained by the previous cell.
            if (output.length !== 0 && output[size - 1].contains(id)) {
                return;
            }
            // Discard any previous cells contained by this cell.
            while (output.length !== 0 && id.contains(output[output.length - 1])) {
                output.splice(output.length - 1, 1);
                // output.remove(output.size() - 1);
            }
            // Check whether the last 3 elements of "output" plus "id" can be
            // collapsed into a single parent cell.
            while (output.length >= 3) {
                size = output.length;
                // A necessary (but not sufficient) condition is that the XOR of the
                // four cells must be zero. This is also very fast to test.
                if ((output[size - 3].id.and(output[size - 2].id).and(output[size - 1].id)).notEquals(id.id)) {
                    break;
                }
                // Now we do a slightly more expensive but exact test. First, compute a
                // mask that blocks out the two bits that encode the child position of
                // "id" with respect to its parent, then check that the other three
                // children all agree with "mask.
                var mask = id.lowestOnBit().shiftLeft(1);
                mask = mask.add(mask.shiftLeft(1)).not();
                // mask = ~(mask + (mask << 1));
                var idMasked = id.id.and(mask);
                if ((output[size - 3].id.and(mask)).notEquals(idMasked)
                    || (output[size - 2].id.and(mask)).notEquals(idMasked)
                    || (output[size - 1].id.and(mask)).notEquals(idMasked) || id.isFace()) {
                    break;
                }
                // Replace four children by their parent cell.
                output.splice(size - 3);
                // output.remove(size - 1);
                // output.remove(size - 2);
                // output.remove(size - 3);
                id = id.parent();
            }
            output.push(id);
        });
        if (output.length < this.size()) {
            this.initRawSwap(output);
            return true;
        }
        return false;
    };
    return S2CellUnion;
}());

/*
 * Copyright 2005 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/**
 * An S2RegionCoverer is a class that allows arbitrary regions to be
 * approximated as unions of cells (S2CellUnion). This is useful for
 * implementing various sorts of search and precomputation operations.
 *
 * Typical usage: {@code S2RegionCoverer coverer; coverer.setMaxCells(5); S2Cap
 * cap = S2Cap.fromAxisAngle(...); S2CellUnion covering;
 * coverer.getCovering(cap, covering); * }
 *
 * This yields a cell union of at most 5 cells that is guaranteed to cover the
 * given cap (a disc-shaped region on the sphere).
 *
 *  The approximation algorithm is not optimal but does a pretty good job in
 * practice. The output does not always use the maximum number of cells allowed,
 * both because this would not always yield a better approximation, and because
 * max_cells() is a limit on how much work is done exploring the possible
 * covering as well as a limit on the final output size.
 *
 *  One can also generate interior coverings, which are sets of cells which are
 * entirely contained within a region. Interior coverings can be empty, even for
 * non-empty regions, if there are no cells that satisfy the provided
 * constraints and are contained by the region. Note that for performance
 * reasons, it is wise to specify a max_level when computing interior coverings
 * - otherwise for regions with small or zero area, the algorithm may spend a
 * lot of time subdividing cells all the way to leaf level to try to find
 * contained cells.
 *
 *  This class is thread-unsafe. Simultaneous calls to any of the getCovering
 * methods will conflict and produce unpredictable results.
 *
 */
var S2RegionCoverer = /** @class */ (function () {
    /**
     * Default constructor, sets all fields to default values.
     */
    function S2RegionCoverer() {
        this.minLevel = 0;
        this.maxLevel = S2CellId.MAX_LEVEL;
        this.levelMod = 1;
        this.maxCells = S2RegionCoverer.DEFAULT_MAX_CELLS;
        this.region = null;
        this.result = [];
        this.candidateQueue = new PriorityQueue();
    }
    // Set the minimum and maximum cell level to be used. The default is to use
    // all cell levels. Requires: max_level() >= min_level().
    //
    // To find the cell level corresponding to a given physical distance, use
    // the S2Cell metrics defined in s2.h. For example, to find the cell
    // level that corresponds to an average edge length of 10km, use:
    //
    // int level = S2::kAvgEdge.GetClosestLevel(
    // geostore::S2Earth::KmToRadians(length_km));
    //
    // Note: min_level() takes priority over max_cells(), i.e. cells below the
    // given level will never be used even if this causes a large number of
    // cells to be returned.
    /**
     * Sets the minimum level to be used.
     */
    S2RegionCoverer.prototype.setMinLevel = function (minLevel) {
        // assert (minLevel >= 0 && minLevel <= S2CellId.MAX_LEVEL);
        this.minLevel = Math.max(0, Math.min(S2CellId.MAX_LEVEL, minLevel));
        return this;
    };
    /**
     * Sets the maximum level to be used.
     */
    S2RegionCoverer.prototype.setMaxLevel = function (maxLevel) {
        // assert (maxLevel >= 0 && maxLevel <= S2CellId.MAX_LEVEL);
        this.maxLevel = Math.max(0, Math.min(S2CellId.MAX_LEVEL, maxLevel));
        return this;
    };
    /**
     * If specified, then only cells where (level - min_level) is a multiple of
     * "level_mod" will be used (default 1). This effectively allows the branching
     * factor of the S2CellId hierarchy to be increased. Currently the only
     * parameter values allowed are 1, 2, or 3, corresponding to branching factors
     * of 4, 16, and 64 respectively.
     */
    S2RegionCoverer.prototype.setLevelMod = function (levelMod) {
        // assert (levelMod >= 1 && levelMod <= 3);
        this.levelMod = Math.max(1, Math.min(3, levelMod));
        return this;
    };
    /**
     * Sets the maximum desired number of cells in the approximation (defaults to
     * kDefaultMaxCells). Note the following:
     *
     * <ul>
     * <li>For any setting of max_cells(), up to 6 cells may be returned if that
     * is the minimum number of cells required (e.g. if the region intersects all
     * six face cells). Up to 3 cells may be returned even for very tiny convex
     * regions if they happen to be located at the intersection of three cube
     * faces.
     *
     * <li>For any setting of max_cells(), an arbitrary number of cells may be
     * returned if min_level() is too high for the region being approximated.
     *
     * <li>If max_cells() is less than 4, the area of the covering may be
     * arbitrarily large compared to the area of the original region even if the
     * region is convex (e.g. an S2Cap or S2LatLngRect).
     * </ul>
     *
     * Accuracy is measured by dividing the area of the covering by the area of
     * the original region. The following table shows the median and worst case
     * values for this area ratio on a test case consisting of 100,000 spherical
     * caps of random size (generated using s2regioncoverer_unittest):
     *
     * <pre>
     * max_cells: 3 4 5 6 8 12 20 100 1000
     * median ratio: 5.33 3.32 2.73 2.34 1.98 1.66 1.42 1.11 1.01
     * worst case: 215518 14.41 9.72 5.26 3.91 2.75 1.92 1.20 1.02
     * </pre>
     */
    S2RegionCoverer.prototype.setMaxCells = function (maxCells) {
        this.maxCells = maxCells;
        return this;
    };
    /**
     * Computes a list of cell ids that covers the given region and satisfies the
     * various restrictions specified above.
     *
     * @param region The region to cover
     * @param covering The list filled in by this method
     */
    S2RegionCoverer.prototype.getCoveringCells = function (region) {
        // Rather than just returning the raw list of cell ids generated by
        // GetCoveringInternal(), we construct a cell union and then denormalize it.
        // This has the effect of replacing four child cells with their parent
        // whenever this does not violate the covering parameters specified
        // (min_level, level_mod, etc). This strategy significantly reduces the
        // number of cells returned in many cases, and it is cheap compared to
        // computing the covering in the first place.
        var tmp = this.getCoveringUnion(region);
        return tmp.denormalize(this.minLevel, this.levelMod);
    };
    /**
     * Computes a list of cell ids that is contained within the given region and
     * satisfies the various restrictions specified above.
     *
     * @param region The region to fill
     * @param interior The list filled in by this method
     */
    S2RegionCoverer.prototype.getInteriorCoveringCells = function (region) {
        var tmp = this.getInteriorCoveringUnion(region);
        return tmp.denormalize(this.minLevel, this.levelMod);
    };
    /**
     * Return a normalized cell union that covers the given region and satisfies
     * the restrictions *EXCEPT* for min_level() and level_mod(). These criteria
     * cannot be satisfied using a cell union because cell unions are
     * automatically normalized by replacing four child cells with their parent
     * whenever possible. (Note that the list of cell ids passed to the cell union
     * constructor does in fact satisfy all the given restrictions.)
     */
    S2RegionCoverer.prototype.getCoveringUnion = function (region, covering) {
        if (covering === void 0) { covering = new S2CellUnion(); }
        this.interiorCovering = false;
        this.getCoveringInternal(region);
        covering.initSwap(this.result);
        return covering;
    };
    /**
     * Return a normalized cell union that is contained within the given region
     * and satisfies the restrictions *EXCEPT* for min_level() and level_mod().
     */
    S2RegionCoverer.prototype.getInteriorCoveringUnion = function (region, covering) {
        if (covering === void 0) { covering = new S2CellUnion(); }
        this.interiorCovering = true;
        this.getCoveringInternal(region);
        covering.initSwap(this.result);
        return covering;
    };
    // /**
    //  * Given a connected region and a starting point, return a set of cells at the
    //  * given level that cover the region.
    //  */
    // public static getSimpleCovering(
    //     region:S2Region , start:S2Point , level:number):S2CellId[] {
    //   S2RegionCoverer.floodFill(region, S2CellId.fromPoint(start).parentL(level));
    // }
    /**
     * If the cell intersects the given region, return a new candidate with no
     * children, otherwise return null. Also marks the candidate as "terminal" if
     * it should not be expanded further.
     */
    S2RegionCoverer.prototype.newCandidate = function (cell) {
        if (!this.region.mayIntersectC(cell)) {
            // console.log("NOT INTERSECTING",this.region);
            return null;
        }
        var isTerminal = false;
        if (cell.level >= this.minLevel) {
            if (this.interiorCovering) {
                if (this.region.containsC(cell)) {
                    isTerminal = true;
                }
                else if (cell.level + this.levelMod > this.maxLevel) {
                    return null;
                }
            }
            else {
                if (cell.level + this.levelMod > this.maxLevel || this.region.containsC(cell)) {
                    isTerminal = true;
                }
            }
        }
        var candidate = new Candidate();
        candidate.cell = cell;
        candidate.isTerminal = isTerminal;
        candidate.numChildren = 0;
        if (!isTerminal) {
            candidate.children = Array.apply(null, new Array(1 << this.maxChildrenShift()));
            // protonew Candidate[1 << this.maxChildrenShift()];
        }
        this.candidatesCreatedCounter++;
        return candidate;
    };
    /** Return the log base 2 of the maximum number of children of a candidate. */
    S2RegionCoverer.prototype.maxChildrenShift = function () {
        return 2 * this.levelMod;
    };
    /**
     * Process a candidate by either adding it to the result list or expanding its
     * children and inserting it into the priority queue. Passing an argument of
     * NULL does nothing.
     */
    S2RegionCoverer.prototype.addCandidate = function (candidate) {
        if (candidate == null) {
            return;
        }
        if (candidate.isTerminal) {
            this.result.push(candidate.cell.id);
            return;
        }
        // assert (candidate.numChildren == 0);
        // Expand one level at a time until we hit min_level_ to ensure that
        // we don't skip over it.
        var numLevels = (candidate.cell.level < this.minLevel) ? 1 : this.levelMod;
        var numTerminals = this.expandChildren(candidate, candidate.cell, numLevels);
        if (candidate.numChildren == 0) ;
        else if (!this.interiorCovering && numTerminals == 1 << this.maxChildrenShift()
            && candidate.cell.level >= this.minLevel) {
            // Optimization: add the parent cell rather than all of its children.
            // We can't do this for interior coverings, since the children just
            // intersect the region, but may not be contained by it - we need to
            // subdivide them further.
            candidate.isTerminal = true;
            this.addCandidate(candidate);
        }
        else {
            // We negate the priority so that smaller absolute priorities are returned
            // first. The heuristic is designed to refine the largest cells first,
            // since those are where we have the largest potential gain. Among cells
            // at the same level, we prefer the cells with the smallest number of
            // intersecting children. Finally, we prefer cells that have the smallest
            // number of children that cannot be refined any further.
            var priority = -((((candidate.cell.level << this.maxChildrenShift()) + candidate.numChildren)
                << this.maxChildrenShift()) + numTerminals);
            this.candidateQueue.add(new QueueEntry(priority, candidate));
            // logger.info("Push: " + candidate.cell.id() + " (" + priority + ") ");
        }
    };
    /**
     * Populate the children of "candidate" by expanding the given number of
     * levels from the given cell. Returns the number of children that were marked
     * "terminal".
     */
    S2RegionCoverer.prototype.expandChildren = function (candidate, cell, numLevels) {
        numLevels--;
        var childCells = cell.subdivide();
        var numTerminals = 0;
        for (var i = 0; i < 4; ++i) {
            if (numLevels > 0) {
                if (this.region.mayIntersectC(childCells[i])) {
                    numTerminals += this.expandChildren(candidate, childCells[i], numLevels);
                }
                continue;
            }
            var child = this.newCandidate(childCells[i]);
            if (child != null) {
                candidate.children[candidate.numChildren++] = child;
                if (child.isTerminal) {
                    ++numTerminals;
                }
            }
        }
        return numTerminals;
    };
    /** Computes a set of initial candidates that cover the given region. */
    S2RegionCoverer.prototype.getInitialCandidates = function () {
        // Optimization: if at least 4 cells are desired (the normal case),
        // start with a 4-cell covering of the region's bounding cap. This
        // lets us skip quite a few levels of refinement when the region to
        // be covered is relatively small.
        if (this.maxCells >= 4) {
            // Find the maximum level such that the bounding cap contains at most one
            // cell vertex at that level.
            var cap = this.region.getCapBound();
            var level = decimal.Decimal.min(S2Projections.MIN_WIDTH.getMaxLevel(cap.angle().radians.times(2)), decimal.Decimal.min(this.maxLevel, S2CellId.MAX_LEVEL - 1)).toNumber();
            if (this.levelMod > 1 && level > this.minLevel) {
                level -= (level - this.minLevel) % this.levelMod;
            }
            // We don't bother trying to optimize the level == 0 case, since more than
            // four face cells may be required.
            if (level > 0) {
                // Find the leaf cell containing the cap axis, and determine which
                // subcell of the parent cell contains it.
                // ArrayList<S2CellId> base = new ArrayList<>(4);
                var id = S2CellId.fromPoint(cap.axis);
                var base = id.getVertexNeighbors(level);
                for (var i = 0; i < base.length; ++i) {
                    this.addCandidate(this.newCandidate(new S2Cell(base[i])));
                }
                return;
            }
        }
        // Default: start with all six cube faces.
        for (var face = 0; face < 6; ++face) {
            this.addCandidate(this.newCandidate(S2RegionCoverer.FACE_CELLS[face]));
        }
    };
    /** Generates a covering and stores it in result. */
    S2RegionCoverer.prototype.getCoveringInternal = function (region) {
        // Strategy: Start with the 6 faces of the cube. Discard any
        // that do not intersect the shape. Then repeatedly choose the
        // largest cell that intersects the shape and subdivide it.
        //
        // result contains the cells that will be part of the output, while the
        // priority queue contains cells that we may still subdivide further. Cells
        // that are entirely contained within the region are immediately added to
        // the output, while cells that do not intersect the region are immediately
        // discarded.
        // Therefore pq_ only contains cells that partially intersect the region.
        // Candidates are prioritized first according to cell size (larger cells
        // first), then by the number of intersecting children they have (fewest
        // children first), and then by the number of fully contained children
        // (fewest children first).
        if (!(this.candidateQueue.size() == 0 && this.result.length == 0)) {
            throw new Error('preconditions are not satisfied');
        }
        // Preconditions.checkState(this.candidateQueue.isEmpty() && this.result.isEmpty());
        this.region = region;
        this.candidatesCreatedCounter = 0;
        this.getInitialCandidates();
        while (this.candidateQueue.size() !== 0 && (!this.interiorCovering || this.result.length < this.maxCells)) {
            var candidate = this.candidateQueue.poll().candidate;
            // logger.info("Pop: " + candidate.cell.id());
            if (candidate.cell.level < this.minLevel || candidate.numChildren == 1
                || this.result.length + (this.interiorCovering ? 0 : this.candidateQueue.size()) + candidate.numChildren
                    <= this.maxCells) {
                // Expand this candidate into its children.
                for (var i = 0; i < candidate.numChildren; ++i) {
                    this.addCandidate(candidate.children[i]);
                }
            }
            else if (this.interiorCovering) ;
            else {
                candidate.isTerminal = true;
                this.addCandidate(candidate);
            }
        }
        this.candidateQueue.clear();
        this.region = null;
    };
    /**
     * By default, the covering uses at most 8 cells at any level. This gives a
     * reasonable tradeoff between the number of cells used and the accuracy of
     * the approximation (see table below).
     */
    S2RegionCoverer.DEFAULT_MAX_CELLS = 8;
    S2RegionCoverer.FACE_CELLS = [0, 1, 2, 3, 4, 5].map(function (face) { return S2Cell.fromFacePosLevel(face, 0, 0); });
    return S2RegionCoverer;
}());
var Candidate = /** @class */ (function () {
    function Candidate() {
    }
    // elements.
    Candidate.prototype.toString = function () {
        return "isTermina: " + this.isTerminal + " - Cell: " + this.cell.toString();
    };
    return Candidate;
}());
var PriorityQueue = /** @class */ (function () {
    function PriorityQueue() {
        this.clear();
    }
    PriorityQueue.prototype.add = function (item) {
        this.items.push(item);
        this.items.sort(function (a, b) { return a.compare(b); });
    };
    PriorityQueue.prototype.clear = function () {
        this.items = [];
    };
    PriorityQueue.prototype.size = function () {
        return this.items.length;
    };
    PriorityQueue.prototype.poll = function () {
        return this.items.splice(0, 1)[0];
    };
    return PriorityQueue;
}());
var QueueEntry = /** @class */ (function () {
    function QueueEntry(id, candidate) {
        this.id = id;
        this.candidate = candidate;
    }
    QueueEntry.prototype.compare = function (other) {
        return this.id < other.id ? 1 : (this.id > other.id ? -1 : 0);
    };
    return QueueEntry;
}());

var Utils = /** @class */ (function () {
    function Utils() {
    }
    /**
     * Calculates a region covering a circle
     * NOTE: The current implementation uses S2Cap while S2Loop would be better (S2Loop is not implemented yet)
     * @param center
     * @param radiusInKM
     * @param points the number of points to calculate. The higher the better precision
     * @returns {S2Region}
     */
    Utils.calcRegionFromCenterRadius = function (center, radiusInKM, points) {
        if (points === void 0) { points = 16; }
        var pointsAtDistance = center.pointsAtDistance(radiusInKM, points);
        var s2Cap = S2Cap.empty().addPoint(center.toPoint());
        // It would be probably enough to add one of the points/2 pair of opposite points in the circle such
        // as (0, points/2). but since this is just a temporary solution lets stick with this as it
        // will come handy when implementing S2Loop.
        pointsAtDistance
            .map(function (p) { return p.toPoint(); })
            .forEach(function (p) {
            s2Cap = s2Cap.addPoint(p);
        });
        return s2Cap;
    };
    return Utils;
}());

exports.Utils = Utils;
exports.Interval = Interval;
exports.MutableInteger = MutableInteger;
exports.R1Interval = R1Interval;
exports.R2Vector = R2Vector;
exports.S1Angle = S1Angle;
exports.S1Interval = S1Interval;
exports.S2 = S2;
exports.S2Metric = S2Metric$$1;
exports.S2Cap = S2Cap;
exports.S2Cell = S2Cell;
exports.S2CellId = S2CellId;
exports.S2CellUnion = S2CellUnion;
exports.S2LatLng = S2LatLng;
exports.S2LatLngRect = S2LatLngRect;
exports.S2Point = S2Point;
exports.S2Projections = S2Projections;
exports.S2RegionCoverer = S2RegionCoverer;
//# sourceMappingURL=export.cjs.js.map
