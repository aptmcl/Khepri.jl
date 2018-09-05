

default_beam_family(
    beam_family_instance(
        backend_family(
            revit => RevitFamily(
                "C:\\ProgramData\\Autodesk\\RVT 2017\\Libraries\\US Metric\\Structural Framing\\Wood\\M_Timber.rfa",
                Dict(:width => "b", :height => "d")),
            autocad => )),
        width=1,
        height=2)
