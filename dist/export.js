import { S2Cap } from "./S2Cap";
export * from './Interval';
export * from './MutableInteger';
export * from './R1Interval';
export * from './R2Vector';
export * from './S1Angle';
export * from './S1Interval';
export * from './S2';
export * from './S2Cap';
export * from './S2Cell';
export * from './S2CellId';
export * from './S2CellUnion';
// export * from './S2EdgeIndex';
// export * from './S2EdgeUtil';
export * from './S2LatLng';
export * from './S2LatLngRect';
// export * from './S2Loop';
export * from './S2Point';
export * from './S2Projections';
export * from './S2RegionCoverer';
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
export { Utils };
//# sourceMappingURL=export.js.map