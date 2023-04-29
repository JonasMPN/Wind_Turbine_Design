function FAST_object = controls(FAST_object, t_max, t_demanded, OMG, A, B, B2, C)
    FAST_object.Control.Torque.Limit = t_max;
    FAST_object.Control.Torque.Demanded = t_demanded;
    FAST_object.Control.Torque.OptGain = OMG;
    FAST_object.Control.Torque.SpeedA = A;
    FAST_object.Control.Torque.SpeedB = B;
    FAST_object.Control.Torque.SpeedB2 = B2;
    FAST_object.Control.Torque.SpeedC = C;
end

