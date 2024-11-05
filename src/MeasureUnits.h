#ifndef MEASUREUNITS_H
#define MEASUREUNITS_H

class MeasureUnits {
public:
    // conversion units /////////////////////////////////////////////////////////////////////////
    // fluid and granular matter are solved in different measure units
    // for a reference, check Feng, Han, Owen, 2007
    // unit length
    double Length;
    //unit time
    double Time;
    //unit density
    double Density;
    // secondary units: speed, acceleration, mass, force, flow rate
    double Volume, Speed, Accel, AngVel, Force, Torque, Mass, KinVisc, DynVisc, Stress, Pressure, FlowRate, Energy;
    // inverted unit length for efficiency
    double invLength;
    // constructor
    constexpr MeasureUnits()
        : Length(1.0)
        , Time(1.0)
        , Density(1.0)
        , Volume(1.0)
        , Speed(1.0)
        , Accel(1.0)
        , AngVel(1.0)
        , Force(1.0)
        , Torque(1.0)
        , Mass(1.0)
        , KinVisc(1.0)
        , DynVisc(1.0)
        , Stress(1.0)
        , Pressure(1.0)
        , FlowRate(0)
        , Energy(0)
        , invLength(0) { }

    void setComposite();
};

#endif  // MEASUREUNITS_H
