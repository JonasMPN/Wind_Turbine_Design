function FAST_object = certification(FAST_object, certification_data)
    FAST_object.CertificationSettings.Wind.Ly = 2*certification_data.new_radius;
    FAST_object.CertificationSettings.Wind.Lz = 2*certification_data.new_radius;
    FAST_object.CertificationSettings.Wind.Class = [certification_data.IEA_wind_class, certification_data.IEA_turbulence_class]';
end

