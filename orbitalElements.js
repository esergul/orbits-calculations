const EARTH_FLATTENING = 1 / 298.257223563; 
const RO = 6378137.0;  // Equatorial radius of the Earth in meters
const RPOL = RO * (1 - EARTH_FLATTENING); // Polar radius of the Earth in meters
const MU = 398600.4418e9; // Earth gravitational parameter

class CartesianVector {
    x;
    y;
    z;
    unit;
    constructor(x, y, z, unit = "m") {
        this.x = x;
        this.y = y;
        this.z = z;
        this.unit = unit;
    }

    norm() {
        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z); 
    }

    dot(rhsVector) {
        return (this.x * rhsVector.x + this.y * rhsVector.y + this.z * rhsVector.z);
    }

    cross(rhsVector) {
        return new CartesianVector(
            this.y * rhsVector.z - this.z * rhsVector.y,
            this.z * rhsVector.x - this.x * rhsVector.z,
            this.x * rhsVector.y - this.y * rhsVector.x
        )
    }
}

class GeodeticElements {
    latitude;
    longitude;
    altitude;
    distanceUnit;
    angleUnit;

    constructor(latitude, longitude, altitude, distanceUnit = "m", angleUnit = "deg") {
        this.latitude = latitude;
        this.longitude = longitude;
        this.altitude = altitude;
        this.distanceUnit = distanceUnit;
        this.angleUnit = angleUnit;
    }
}

class KeplerianElements {
    semiMajorAxis;
    eccentricity;
    inclination;
    rightAscensionOfTheAscendingNode;
    argumentOfPerigee;
    trueAnomaly;
    distanceUnit;
    angleUnit;

    constructor(a, e, i, raan, w, ta, distanceUnit = "m", angleUnit = "deg") {
        this.semiMajorAxis = a;
        this.eccentricity = e;
        this.inclination = i;
        this.rightAscensionOfTheAscendingNode = raan;
        this.argumentOfPerigee = w;
        this.trueAnomaly = ta;
        this.distanceUnit = distanceUnit;
        this.angleUnit = angleUnit;
    }

    compare(rhs) {
        if (this.angleUnit != rhs.angleUnit || this.distanceUnit != rhs.distanceUnit) {
            console.log("Units mismatch, cant compare");
            return -1;
        }
        const err = 1e-4;
        return Math.abs(this.semiMajorAxis - rhs.semiMajorAxis) < err &&
        Math.abs(this.eccentricity - rhs.eccentricity) < err &&
        Math.abs(this.inclination - rhs.inclination) < err &&
        Math.abs(this.rightAscensionOfTheAscendingNode - rhs.rightAscensionOfTheAscendingNode) < err &&
        Math.abs(this.argumentOfPerigee - rhs.argumentOfPerigee) < err &&
        Math.abs(this.trueAnomaly - rhs.trueAnomaly) < err;
    }
}

function cartesianToKeplerianElements(R, V) {
    const EPS = 1e-10; 
    const PI = Math.PI;

    const r = R.norm();
    const v = V.norm();
    const vr = R.dot(V) / r;

    const H = R.cross(V);
    const h = H.norm();

    const incl = Math.acos(H.z / h);

    const N = new CartesianVector(0, 0, 1).cross(H);
    const n = N.norm();

    let RA = 0;
    if (n !== 0) {
        RA = Math.acos(N.x / n);
        if (N.y < 0)
            RA = 2 * PI - RA;
    }

    const muReciprocal = 1 / MU;
    const A = (v * v - MU / r);
    const E = new CartesianVector(
        muReciprocal * (A * R.x - r * vr * V.x),
        muReciprocal * (A * R.y - r * vr * V.y),
        muReciprocal * (A * R.z - r * vr * V.z)
    );
    const e = E.norm();

    let w = 0;
    if (n !== 0 && e > EPS) {
        w = Math.acos((N.dot(E) / (n * e)));
        if (E.z < 0)
            w = 2 * PI - w;
    }

    let TA;
    let cp;

    if (e > EPS) {
        TA = Math.acos(E.dot(R) / (e * r));
        if (vr < 0)
            TA = 2 * PI - TA;
    } else {
        cp = N.cross(R);
        TA = (cp.z) >= 0 ? acos(N.dot(R) / (n *r)) : 2*PI - acos(N.dot(R) / (n *r));    
    }

    if (isNaN(TA))
        TA = 0;

    const a = h * h / MU / (1 - e * e);

    return new KeplerianElements (a, e, incl * 180 / PI, w * 180 / PI, RA * 180 / PI, TA * 180 / PI);
}

function cartesianToGeodeticPosition(R) {
    const e_sq = EARTH_FLATTENING * (2 - EARTH_FLATTENING);
    const eps = e_sq / (1.0 - e_sq);
    const p = Math.sqrt(R.x * R.x + R.y * R.y);
    const q = Math.atan2(R.z * RO, p * RPOL);
    const sin_q_3 = Math.pow(Math.sin(q), 3);
    const cos_q_3 = Math.pow(Math.cos(q), 3);
    const phi = Math.atan2(R.z + eps * RPOL * sin_q_3, p - e_sq * RO * cos_q_3);
    const lambda = Math.atan2(R.y, R.x);
    const v = RO / Math.sqrt(1 - e_sq * Math.sin(phi) * Math.sin(phi));
    const h = p / Math.cos(phi) - v;
  
    return new GeodeticElements(phi * 180 / Math.PI, lambda * 180 / Math.PI, h);
}

function test() {
    const cartesianEciPositionVector = new CartesianVector(-5018330.730527248, -4842014.696922994, 378279.861098418);
    const cartesianEciVelocityVector =  new CartesianVector(-387.4431638836936, 990.9497342536506, 7477.068695495851);
    const cartesianEcefPositionVector = new CartesianVector(-6082840.226777075, 3410994.051399341, 368808.42765045795);
    const geodeticPosition = new GeodeticElements(3.045785642831048, 150.71817446508888, 605607.1262781937);
    const keplerianElements = new KeplerianElements(6979122.709428061, 8.124206313827294E-4, 97.4611130476609, 
        146.72136147502, 224.38262613453847, 216.41019492014885);

    let keplerianConversionResult = cartesianToKeplerianElements(cartesianEciPositionVector, cartesianEciVelocityVector);
    if (keplerianConversionResult.compare(keplerianElements)) {
        console.log("Keplerian conversion OK");
    } else {
        console.log("Keplerian conversion error");
    }

    let geodeticConversionResult = cartesianToGeodeticPosition(cartesianEcefPositionVector);
    if (geodeticPosition.latitude.toPrecision(3) == geodeticConversionResult.latitude.toPrecision(3) &&
        geodeticPosition.longitude.toPrecision(3) == geodeticConversionResult.longitude.toPrecision(3) &&
        geodeticPosition.altitude.toPrecision(3) == geodeticConversionResult.altitude.toPrecision(3)) {
        console.log("Geodetic conversion OK");
    } else {
        console.log("Geodetic conversion error");
    }
}

test();