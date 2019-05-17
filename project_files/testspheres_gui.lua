mmCreateView("testspheres", "GUIView", "::gui")
mmCreateModule("View3D", "::testview")
mmCreateModule("SimpleSphereRenderer", "::rnd")
mmSetParamValue("rnd::renderMode", "Simple")
mmCreateModule("TestSpheresDataSource", "::dat")
mmCreateCall("CallRenderView", "::gui::renderview", "::testview::render")
mmCreateCall("CallRender3D", "::testview::rendering", "::rnd::rendering")
mmCreateCall("MultiParticleDataCall", "::rnd::getData", "::dat::getData")
