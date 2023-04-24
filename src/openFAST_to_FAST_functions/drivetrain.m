function FAST_object = drivetrain(FAST_object, generator_efficiency, gearbox_ratio)
    FAST_object.Drivetrain.Generator.Efficiency = generator_efficiency;
    FAST_object.Drivetrain.Gearbox.Ratio = gearbox_ratio;
end

