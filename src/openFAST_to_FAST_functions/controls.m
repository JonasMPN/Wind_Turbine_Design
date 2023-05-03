function FAST_object = controls(FAST_object, control_data)
    FAST_object.Control.Torque.Demanded = control_data.torque_demanded;
    FAST_object.Control.Torque.Min = control_data.torque_min;
    FAST_object.Control.Torque.Slewrate = control_data.torque_demanded/2;
    FAST_object.Control.Torque.Limit = control_data.torque_max;
    FAST_object.Control.Torque.OptGain = control_data.opt_mode_gain;
    FAST_object.Control.Torque.SpeedA = control_data.omega_A;
    FAST_object.Control.Torque.SpeedB = control_data.omega_B;
    FAST_object.Control.Torque.SpeedB2 = control_data.omega_B2;
    FAST_object.Control.Torque.SpeedC = control_data.omega_C;
end

