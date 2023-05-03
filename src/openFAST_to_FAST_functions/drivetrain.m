function FAST_object = drivetrain(FAST_object, drivetrain_data)
    FAST_object.Drivetrain.Generator.Efficiency = drivetrain_data.gen_efficiency;
    FAST_object.Drivetrain.Gearbox.Ratio = drivetrain_data.gearbox_ratio;
    FAST_object.Drivetrain.Generator.HSSInertia = drivetrain_data.gen_inertia;
end

