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

    void setComposite() {
        Volume = Length * Length * Length;
        Speed = Length / Time; //+2
        Accel = Length / Time / Time; // +6
        AngVel = 1.0 / Time;
        KinVisc = Length * Length / Time; // 0
        DynVisc = Density * Length * Length / Time; // 0
        Force = Density * Length * Length * Length * Length / Time / Time; // +3
        Torque = Density * Length * Length * Length * Length * Length / Time / Time; // +1
        Mass = Density * Length * Length * Length; // -3
        Stress = Density * Length * Length / Time / Time;
        Pressure = Density * Length * Length / Time / Time;
        FlowRate = Density * Length * Length * Length / Time;
        Energy = Density * Length * Length * Length * Length * Length / Time / Time;
        invLength = 1.0 / Length;
    }
};

#endif  // MEASUREUNITS_H
