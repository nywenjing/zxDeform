#include "zxsimviewer.h"
#include "zxmath.h"
#include "zxneohookeanmaterial.h"
#include "zxnonlinearfem_forcemodel_sparse.h"
#include "zxtimestepperbackwardeuler.h"
#include "zxcontactconstraint.h"
#include "igl/per_vertex_normals.h"
#include <GL/gl.h>

zxSimViewer::zxSimViewer()
{

    zxContactPoint::Ptr cp = zxContactPoint::create();
    zxTestConstraint* c = new zxTestConstraint(cp);

    m_world = zxSimWorld::create();

    init_femsim0();


}

void zxSimViewer::init()
{

}

void check_validility(zxBVHNode::Ptr node)
{
    if(node->get_aabb() == nullptr)
        std::cout<<"error"<<std::endl;

    if(node->is_leaf())
        return;
    check_validility(node->get_left());
    check_validility(node->get_right());

}

void zxSimViewer::init_femsim0()
{

    std::string volfilename = "../zxDeform/data/torus.1";
    zxSolidMesh::Ptr vmesh = zxTetrahedralMesh::create(volfilename);
    zxRenderMesh::Ptr rmesh = zxRenderMesh::create(vmesh);
    zxCollisionMesh::Ptr cmesh  = rmesh;
    zxMaterial::Ptr material = zxNeoHookeanMaterial::create(1e4,0.0,1e3);
    zxNonlinearFEM_ForceModel_Sparse::Ptr forcemodel = zxNonlinearFEM_ForceModel_Sparse::create(vmesh,material);
    zxTimeStepperBackwardEuler::Ptr stepper = zxTimeStepperBackwardEuler::create(forcemodel);
    zxBody::Ptr body = zxBody::create();
    body->set_render_mesh(rmesh);
    body->set_collision_mesh(cmesh);
    body->set_stepper(stepper);
    m_world->add_body(body);

    body->get_bvh_tree()->refit();
    check_validility(body->get_bvh_tree()->get_root());

    std::string planeName = "../zxDeform/data/plane.obj";
    zxRenderMesh::Ptr pmesh = zxRenderMesh::create(planeName);
    zxBody::Ptr plane = zxBody::create();
    plane->set_render_mesh(pmesh);
    plane->set_collision_mesh(pmesh);

    m_world->add_body(plane);

    m_rad = 6;

}

void zxSimViewer::render()
{
    //m_rad = 6;
    mat4d mv_mat = get_modelview_matrix();//zx_gl_lookAt(m_cen - m_forward * m_rad * 2,m_cen,m_up);
    mat4d proj_mat = get_perspective_matrix();//zx_gl_perspective(60.0,width() * 1.0/ height(),0.1,300);

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixd(mv_mat.data());
    glMatrixMode(GL_PROJECTION);
    glLoadMatrixd(proj_mat.data());

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glShadeModel (GL_SMOOTH);
    glClearColor (1.0, 1.0, 1.0, 1.0);

    for(size_t ib = 0; ib < m_world->get_num_bodies(); ib++)
    {
        zxBody::Ptr body = m_world->get_body(ib);
        zxRenderMesh::Ptr r_mesh = body->get_render_mesh();

        Eigen::MatrixXd V,N;
        Eigen::MatrixXi F;
        r_mesh->get_Matrix_Format(V,F);
        igl::per_vertex_normals(V,F,N);

        glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,r_mesh->m_ambient.data());
        glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,r_mesh->m_diffuse.data());
        glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,r_mesh->m_specular.data());
        glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,r_mesh->m_shiness);

        glBegin(GL_TRIANGLES);
        for(int el = 0; el < F.rows(); el++)
            for(int i = 0; i < 3; i++)
            {
                vec3d x = V.row(F(el,i));
                vec3d n = N.row(F(el,i));

                glNormal3dv(n.data());
                glVertex3dv(x.data());
            }

        glEnd();
    }

}

void zxSimViewer::animate()
{
    if(m_isAnimated)
        m_world->do_time_step();
}
