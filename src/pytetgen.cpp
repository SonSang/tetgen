#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <string>
#include <cassert>
#include "../tetgen.h"

namespace py = pybind11;


PYBIND11_MODULE(pytetgen, m) {
    m.doc() = "python binding of tetgen";
    
    /*
    tetgenbehavior
    */
    py::class_<tetgenbehavior>(m, "tetgenbehavior")
        .def(py::init<>())
        .def("parse_commandline", [](tetgenbehavior &self, const std::string& commands) {
            self.parse_commandline(const_cast<char*>(commands.c_str()));
        })
        .def_readonly("plc", &tetgenbehavior::plc)
        .def_readonly("psc", &tetgenbehavior::psc)
        .def_readonly("refine", &tetgenbehavior::refine)
        .def_readonly("quality", &tetgenbehavior::quality)
        .def_readonly("nobisect", &tetgenbehavior::nobisect)
        .def_readonly("cdt", &tetgenbehavior::cdt)
        .def_readonly("cdtrefine", &tetgenbehavior::cdtrefine)
        .def_readonly("coarsen", &tetgenbehavior::coarsen)
        .def_readonly("weighted", &tetgenbehavior::weighted)
        .def_readonly("brio_hilbert", &tetgenbehavior::brio_hilbert)
        .def_readonly("flipinsert", &tetgenbehavior::flipinsert)
        .def_readonly("metric", &tetgenbehavior::metric)
        .def_readonly("varvolume", &tetgenbehavior::varvolume)
        .def_readonly("fixedvolume", &tetgenbehavior::fixedvolume)
        .def_readonly("regionattrib", &tetgenbehavior::regionattrib)
        .def_readonly("insertaddpoints", &tetgenbehavior::insertaddpoints)
        .def_readonly("diagnose", &tetgenbehavior::diagnose)
        .def_readonly("convex", &tetgenbehavior::convex)
        .def_readonly("nomergefacet", &tetgenbehavior::nomergefacet)
        .def_readonly("nomergevertex", &tetgenbehavior::nomergevertex)
        .def_readonly("noexact", &tetgenbehavior::noexact)
        .def_readonly("nostaticfilter", &tetgenbehavior::nostaticfilter)
        .def_readonly("zeroindex", &tetgenbehavior::zeroindex)
        .def_readonly("facesout", &tetgenbehavior::facesout)
        .def_readonly("edgesout", &tetgenbehavior::edgesout)
        .def_readonly("neighout", &tetgenbehavior::neighout)
        .def_readonly("voroout", &tetgenbehavior::voroout)
        .def_readonly("meditview", &tetgenbehavior::meditview)
        .def_readonly("vtkview", &tetgenbehavior::vtkview)
        .def_readonly("vtksurfview", &tetgenbehavior::vtksurfview)
        .def_readonly("nobound", &tetgenbehavior::nobound)
        .def_readonly("nonodewritten", &tetgenbehavior::nonodewritten)
        .def_readonly("noelewritten", &tetgenbehavior::noelewritten)
        .def_readonly("nofacewritten", &tetgenbehavior::nofacewritten)
        .def_readonly("noiterationnum", &tetgenbehavior::noiterationnum)
        .def_readonly("nojettison", &tetgenbehavior::nojettison)
        .def_readonly("docheck", &tetgenbehavior::docheck)
        .def_readonly("quiet", &tetgenbehavior::quiet)
        .def_readonly("nowarning", &tetgenbehavior::nowarning)
        .def_readonly("verbose", &tetgenbehavior::verbose)
        .def_readonly("vertexperblock", &tetgenbehavior::vertexperblock)
        .def_readonly("tetrahedraperblock", &tetgenbehavior::tetrahedraperblock)
        .def_readonly("shellfaceperblock", &tetgenbehavior::shellfaceperblock)
        .def_readonly("supsteiner_level", &tetgenbehavior::supsteiner_level)
        .def_readonly("addsteiner_algo", &tetgenbehavior::addsteiner_algo)
        .def_readonly("coarsen_param", &tetgenbehavior::coarsen_param)
        .def_readonly("weighted_param", &tetgenbehavior::weighted_param)
        .def_readonly("fliplinklevel", &tetgenbehavior::fliplinklevel)
        .def_readonly("flipstarsize", &tetgenbehavior::flipstarsize)
        .def_readonly("fliplinklevelinc", &tetgenbehavior::fliplinklevelinc)
        .def_readonly("opt_max_flip_level", &tetgenbehavior::opt_max_flip_level)
        .def_readonly("opt_scheme", &tetgenbehavior::opt_scheme)
        .def_readonly("opt_iterations", &tetgenbehavior::opt_iterations)
        .def_readonly("smooth_cirterion", &tetgenbehavior::smooth_cirterion)
        .def_readonly("smooth_maxiter", &tetgenbehavior::smooth_maxiter)
        .def_readonly("delmaxfliplevel", &tetgenbehavior::delmaxfliplevel)
        .def_readonly("order", &tetgenbehavior::order)
        .def_readonly("reversetetori", &tetgenbehavior::reversetetori)
        .def_readonly("steinerleft", &tetgenbehavior::steinerleft)
        .def_readonly("unflip_queue_limit", &tetgenbehavior::unflip_queue_limit)
        .def_readonly("no_sort", &tetgenbehavior::no_sort)
        .def_readonly("hilbert_order", &tetgenbehavior::hilbert_order)
        .def_readonly("hilbert_limit", &tetgenbehavior::hilbert_limit)
        .def_readonly("brio_threshold", &tetgenbehavior::brio_threshold)
        .def_readonly("brio_ratio", &tetgenbehavior::brio_ratio)
        .def_readonly("epsilon", &tetgenbehavior::epsilon)
        .def_readonly("facet_separate_ang_tol", &tetgenbehavior::facet_separate_ang_tol)
        .def_readonly("collinear_ang_tol", &tetgenbehavior::collinear_ang_tol)
        .def_readonly("facet_small_ang_tol", &tetgenbehavior::facet_small_ang_tol)
        .def_readonly("maxvolume", &tetgenbehavior::maxvolume)
        .def_readonly("maxvolume_length", &tetgenbehavior::maxvolume_length)
        .def_readonly("minratio", &tetgenbehavior::minratio)
        .def_readonly("opt_max_asp_ratio", &tetgenbehavior::opt_max_asp_ratio)
        .def_readonly("opt_max_edge_ratio", &tetgenbehavior::opt_max_edge_ratio)
        .def_readonly("mindihedral", &tetgenbehavior::mindihedral)
        .def_readonly("optmaxdihedral", &tetgenbehavior::optmaxdihedral)
        .def_readonly("metric_scale", &tetgenbehavior::metric_scale)
        .def_readonly("smooth_alpha", &tetgenbehavior::smooth_alpha)
        .def_readonly("coarsen_percent", &tetgenbehavior::coarsen_percent)
        .def_readonly("elem_growth_ratio", &tetgenbehavior::elem_growth_ratio)
        .def_readonly("refine_progress_ratio", &tetgenbehavior::refine_progress_ratio)
        ;

    /*
    tetgenio
    */

    // polygon
    py::class_<tetgenio::polygon>(m, "polygon")
        .def(py::init<>())
        .def("set_vertexlist", [](tetgenio::polygon &self, py::array_t<int> vertexlist) {
            auto r = vertexlist.unchecked<1>();
            self.numberofvertices = r.shape(0);
            self.vertexlist = new int[r.shape(0)];
            for (int i = 0; i < r.shape(0); i++) {
                self.vertexlist[i] = r(i);
            }
        })
        .def_readonly("numberofvertices", &tetgenio::polygon::numberofvertices)
        .def_readonly("vertexlist", &tetgenio::polygon::vertexlist)
        ;

    // facet
    py::class_<tetgenio::facet>(m, "facet")
        .def(py::init<>())
        .def("set_polygonlist", [](tetgenio::facet &self, std::vector<tetgenio::polygon> polygonlist) {
            self.numberofpolygons = polygonlist.size();
            self.polygonlist = new tetgenio::polygon[self.numberofpolygons];
            for (int i = 0; i < self.numberofpolygons; i++) {
                self.polygonlist[i] = polygonlist[i];
            }
        })
        .def_readonly("numberofpolygons", &tetgenio::facet::numberofpolygons)
        .def_readonly("polygonlist", &tetgenio::facet::polygonlist)
        .def_readonly("numberofholes", &tetgenio::facet::numberofholes)
        .def_readonly("holelist", &tetgenio::facet::holelist)
        ;

    // tetgenio
    py::class_<tetgenio>(m, "tetgenio")
        .def(py::init<>())
        .def("set_pointlist", [](tetgenio &self, py::array_t<REAL> pointlist) {
            auto r = pointlist.unchecked<2>();
            self.numberofpoints = r.shape(0);
            self.pointlist = new REAL[r.shape(0) * 3];
            for (int i = 0; i < r.shape(0); i++) {
                for (int j = 0; j < 3; j++) {
                    self.pointlist[i * 3 + j] = r(i, j);
                }
            }
        })
        .def("set_pointmarkerlist", [](tetgenio &self, std::vector<int> pointmarkerlist) {
            self.pointmarkerlist = new int[pointmarkerlist.size()];
            for (int i = 0; i < int(pointmarkerlist.size()); i++) {
                self.pointmarkerlist[i] = pointmarkerlist[i];
            }
        })
        .def("set_facetlist_from_trimesh", [](tetgenio &self, py::array_t<int> facetlist) {
            auto r = facetlist.unchecked<2>();
            assert (r.shape(1) == 3);   // should be triangular mesh;

            self.numberoffacets = r.shape(0);
            self.facetlist = new tetgenio::facet[self.numberoffacets];
            for (int i = 0; i < self.numberoffacets; i++) {
                self.facetlist[i].numberofpolygons = 1;
                self.facetlist[i].polygonlist = new tetgenio::polygon[1];
                self.facetlist[i].polygonlist[0].numberofvertices = 3;
                self.facetlist[i].polygonlist[0].vertexlist = new int[3];
                for (int j = 0; j < 3; j++) {
                    self.facetlist[i].polygonlist[0].vertexlist[j] = r(i, j);
                }
            }
        })
        .def("set_facetlist", [](tetgenio &self, std::vector<tetgenio::facet> facetlist) {
            self.numberoffacets = facetlist.size();
            self.facetlist = new tetgenio::facet[self.numberoffacets];
            for (int i = 0; i < self.numberoffacets; i++) {
                self.facetlist[i] = facetlist[i];
            }
        })
        .def("set_facetmarkerlist", [](tetgenio &self, std::vector<int> facetmarkerlist) {
            self.facetmarkerlist = new int[facetmarkerlist.size()];
            for (int i = 0; i < int(facetmarkerlist.size()); i++) {
                self.facetmarkerlist[i] = facetmarkerlist[i];
            }
        })
        .def("get_pointlist", [](tetgenio &self) {
            py::array_t<REAL> result({self.numberofpoints, 3});
            auto r = result.mutable_unchecked<2>();
            for (int i = 0; i < self.numberofpoints; i++) {
                for (int j = 0; j < 3; j++) {
                    r(i, j) = self.pointlist[i * 3 + j];
                }
            }
            return result;
        })
        .def("get_trifacelist", [](tetgenio &self) {
            py::array_t<int> result({self.numberoftrifaces, 3});
            auto r = result.mutable_unchecked<2>();
            for (int i = 0; i < self.numberoftrifaces; i++) {
                for (int j = 0; j < 3; j++) {
                    r(i, j) = self.trifacelist[i * 3 + j];
                }
            }
            return result;
        })
        .def("get_trifacemarkerlist", [](tetgenio &self) {
            py::array_t<int> result({self.numberoftrifaces});
            auto r = result.mutable_unchecked<1>();
            for (int i = 0; i < self.numberoftrifaces; i++) {
                r(i) = self.trifacemarkerlist[i];
            }
            return result;
        })
        .def("get_tetrahedronlist", [](tetgenio &self) {
            py::array_t<int> result({self.numberoftetrahedra, 4});
            auto r = result.mutable_unchecked<2>();
            for (int i = 0; i < self.numberoftetrahedra; i++) {
                for (int j = 0; j < 4; j++) {
                    r(i, j) = self.tetrahedronlist[i * 4 + j];
                }
            }
            return result;
        })
        .def_readwrite("firstnumber", &tetgenio::firstnumber)
        
        .def_readonly("pointlist", &tetgenio::pointlist)
        .def_readonly("pointmarkerlist", &tetgenio::pointmarkerlist)
        .def_readonly("numberofpoints", &tetgenio::numberofpoints)
        
        .def_readonly("tetrahedronlist", &tetgenio::tetrahedronlist)
        .def_readonly("numberoftetrahedra", &tetgenio::numberoftetrahedra)
        
        .def_readonly("facetlist", &tetgenio::facetlist)
        .def_readonly("facetmarkerlist", &tetgenio::facetmarkerlist)
        .def_readonly("numberoffacets", &tetgenio::numberoffacets)

        .def_readonly("trifacelist", &tetgenio::trifacelist)
        .def_readonly("trifacemarkerlist", &tetgenio::trifacemarkerlist)
        .def_readonly("numberoftrifaces", &tetgenio::numberoftrifaces)

        .def_readonly("edgelist", &tetgenio::edgelist)
        .def_readonly("edgemarkerlist", &tetgenio::edgemarkerlist)
        .def_readonly("numberofedges", &tetgenio::numberofedges)
        ;

    /*
    tetrahedralize
    */
    m.def("tetrahedralize", [](tetgenbehavior &b, tetgenio &in, tetgenio &out) {
        tetrahedralize(&b, &in, &out);
    });
}