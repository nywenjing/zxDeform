#include "zxalmsimulator.h"

zxALMSimulator::zxALMSimulator()
{
    m_residual_tol = 0.01;
    m_disp_tol = 1e-6;
    m_energy_tol = 0.01;
}

void zxALMSimulator::add_forcemodel(zxForceModel::Ptr forcemodel)
{
    m_forceModels.push_back(forcemodel);

    zxMesh::Ptr mesh = forcemodel->get_mesh();
    std::vector<zxNode::Ptr>& nodes = mesh->get_nodes();
    m_nodes.insert(m_nodes.end(),nodes.begin(),nodes.end());
}

void zxALMSimulator::add_contactInterface(zxALMContactInterface::Ptr contactInterface)
{
    m_contactInterfaces.push_back(contactInterface);
}

void zxALMSimulator::init_simulator()
{


}

void zxALMSimulator::init_step()
{
    m_numDofs = 0;

    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        zxNode::Ptr node = m_nodes[i];

        for(size_t j = 0; j < 3; j++)
        {
            if(node->m_bc[j] == zxFree)
            {
                node->m_dof_id[j] = m_numDofs++;
            }
            else if(node->m_bc[j] == zxFixed)
                node->m_dof_id[j] = -1;
            else if(node->m_bc[j] == zxPrescribed)
                node->m_dof_id[j] = -2 - m_numDofs++;
        }
    }

    m_ui = m_Ui = m_Ut = m_R0 = m_R1 = Eigen::VectorXd(m_numDofs);
    m_ui.setZero(); m_Ui.setZero(); m_Ut.setZero();

    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        zxNode* node = m_nodes[i].get();

        node->rp = node->rt;

        for(size_t j = 0; j < 3; j++)
        {
            int eid = node->m_dof_id[j];
            if(eid >= 0)
                m_Ut[eid] = node->rt[j] - node->r0[j];
        }
    }
    m_stiffness_matrix.resize(m_numDofs,m_numDofs);

    m_is_prescribed_dofs.resize(m_numDofs);
    m_is_prescribed_dofs.assign(m_is_prescribed_dofs.size(),false);

    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        zxNode::Ptr node = m_nodes[i];

        for(size_t j = 0; j < 3; j++)
        {
            if(node->m_bc[j] == zxPrescribed)
            {
                int id = -node->m_dof_id[j] - 2;
                m_is_prescribed_dofs[id] = true;
            }
        }
    }

    update_contact();
    update_stress();
}

void zxALMSimulator::do_simulate()
{
    if(m_numDofs == 0)
        return;

    init_step();
    compute_stiffness_matrix();
    compute_residual(m_R0);

    bool converge  = false;

    real rnorm,enorm,unorm,Unorm,initRnorm,initEnorm,initUnorm;

    size_t niter = 0;

    m_R0 -= m_Fd;
    while(!converge)
    {
        converge = true;
        m_linearSolver.compute(m_stiffness_matrix);
        m_ui = m_linearSolver.solve(m_R0);

        if(niter == 0)
        {
            initRnorm = m_R0.squaredNorm();
            initEnorm = m_ui.dot(m_R0);
            initUnorm = m_ui.squaredNorm();
        }

        real s = do_line_search(1.0);
        m_ui *= s;
        m_Ui += m_ui;

        rnorm = m_R1.squaredNorm();
        unorm = m_ui.squaredNorm();
        Unorm = m_Ui.squaredNorm();
        enorm = std::abs(m_ui.dot(m_R1));

        if(rnorm > m_residual_tol * initRnorm) converge = false;
        if(unorm > m_disp_tol * initUnorm) converge = false;
        if(enorm > m_energy_tol * initEnorm) converge = false;

        printf("\tstep from line search         = %lf\n", s);

        printf("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
        printf("\t   residual         %15le %15le %15le \n", initRnorm, rnorm, m_residual_tol*initRnorm);
        printf("\t   energy           %15le %15le %15le \n", initEnorm, enorm, m_energy_tol*initEnorm);
        printf("\t   displacement     %15le %15le %15le \n", initUnorm, unorm ,m_disp_tol*initUnorm );

        fflush(stdout);

        if(!converge)
        {
            m_ui.setZero();
            compute_stiffness_matrix();
            m_R0 = m_R1;
        }
        else
        {
            converge = alm_augment();

            if(!converge)
            {
                update_stress();
                compute_residual(m_R0);
                compute_stiffness_matrix();
            }
        }

        niter++;
    };
}

void zxALMSimulator::update_deformation(Eigen::VectorXd &u)
{
    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        zxNode* node = m_nodes[i].get();

        for(size_t j = 0; j < 3; j++)
        {
            int eid = node->m_dof_id[j];
            if(eid >= 0)
                node->rt[j] = node->r0[j] + m_Ut[eid] + m_Ui[eid] + u[eid];

            if(-eid - 2 >= 0)
            {
                eid = -eid - 2;
                node->rt[j] = node->rl[j];
            }
        }
    }

}

void zxALMSimulator::update_stress()
{
    for(size_t ib = 0; ib < m_forceModels.size(); ib++)
    {
        zxForceModel::Ptr forcemodel = m_forceModels[ib];
        zxMesh::Ptr mesh = forcemodel->get_mesh();
        Eigen::VectorXd u(mesh->get_num_nodes() * 3);

        for(size_t i = 0; i < mesh->get_num_nodes(); i++)
            for(size_t j = 0; j < 3; j++)
                u[3 * i + j] = mesh->get_node(i)->rt[j] - mesh->get_node(i)->r0[j];
        forcemodel->updatePosition(u);
    }

}

void zxALMSimulator::update_contact()
{
    for(size_t ic = 0; ic < m_contactInterfaces.size(); ic++)
    {
        zxALMContactInterface* cinterface = m_contactInterfaces[ic].get();
        cinterface->update_contact();
    }

}

void zxALMSimulator::compute_residual(Eigen::VectorXd &R)
{
    R.setZero();

    for(size_t ni = 0; ni < m_nodal_forces.size(); ni++)
    {
        zxNodalForce* nforce = m_nodal_forces[ni].get();
        for(size_t i = 0; i < 3; i++)
            if(nforce->m_node->m_dof_id[i] >= 0)
                R[nforce->m_node->m_dof_id[i]] += nforce->m_force[i];
    }

    for(size_t ib = 0; ib < m_forceModels.size(); ib++)
    {
        zxForceModel::Ptr forcemodel = m_forceModels[ib];
        zxMesh::Ptr mesh = forcemodel->get_mesh();
        Eigen::VectorXd force(mesh->get_num_nodes() * 3);
        force.setZero();
        forcemodel->computeForce(force);

        for(size_t i = 0; i < mesh->get_num_nodes(); i++)
            for(size_t j = 0; j < 3; j++)
                if(mesh->get_node(i)->m_bc[j] == zxFree)
                    R[mesh->get_node(i)->m_dof_id[j]] -= force[3 * i + j];
    }

    for(size_t ic = 0; ic < m_contactInterfaces.size(); ic++)
    {
        zxALMContactInterface* cinterface = m_contactInterfaces[ic].get();
        cinterface->compute_contact_force(R);
    }


}

void zxALMSimulator::reform_stiffness_topology()
{
    std::list<Eigen::Triplet<real>> stiffEntries;
    for(size_t ib = 0; ib < m_forceModels.size(); ib++)
    {
        zxForceModel::Ptr forcemodel = m_forceModels[ib];
        zxMesh::Ptr mesh = forcemodel->get_mesh();
        const zxSparseMatrix&   sK = forcemodel->getStiffMatrixTopology();

        for (int k=0; k<sK.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(sK,k); it; ++it)
            {
                int rid = it.row();
                int cid = it.col();
                zxNode::Ptr node0 = mesh->get_node(rid/3);
                zxNode::Ptr node1 = mesh->get_node(cid/3);

                int dof0 = node0->m_dof_id[rid%3];
                int dof1 = node1->m_dof_id[cid%3];

                if(dof0 < -1)
                    dof0 = - dof0 - 2;
                if(dof1 < -1)
                    dof1 = - dof1 - 2;

                if(dof0 >= 0 && dof1 >= 0)
                    stiffEntries.push_back(Eigen::Triplet<real>(dof0,dof1,0.0));
            }
    }

    for(size_t ic = 0; ic < m_contactInterfaces.size(); ic++)
    {
        zxALMContactInterface* cinterface = m_contactInterfaces[ic].get();

        std::list<Eigen::Triplet<real>> clist = cinterface->get_stiffness_triplets();

        stiffEntries.insert(stiffEntries.end(),clist.begin(),clist.end());
    }

    m_stiffness_matrix.setFromTriplets(stiffEntries.begin(),stiffEntries.end());

}

void zxALMSimulator::compute_stiffness_matrix()
{
    reform_stiffness_topology();
    for(size_t ib = 0; ib < m_forceModels.size(); ib++)
    {
        zxForceModel::Ptr forcemodel = m_forceModels[ib];
        zxMesh::Ptr mesh = forcemodel->get_mesh();
        zxSparseMatrix   sK = forcemodel->getStiffMatrixTopology();
        forcemodel->computeTangent(sK);

        for (int k=0; k<sK.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(sK,k); it; ++it)
            {
                int rid = it.row();
                int cid = it.col();
                zxNode::Ptr node0 = mesh->get_node(rid/3);
                zxNode::Ptr node1 = mesh->get_node(cid/3);

                int dof0 = node0->m_dof_id[rid%3];
                int dof1 = node1->m_dof_id[cid%3];

                if(dof0 < -1)
                    dof0 = - dof0 - 2;
                if(dof1 < -1)
                    dof1 = - dof1 - 2;

                if(dof0 >= 0 && dof1 >= 0)
                    m_stiffness_matrix.coeffRef(dof0,dof1) += it.value();
            }
    }

    for(size_t ic = 0; ic < m_contactInterfaces.size(); ic++)
    {
        zxALMContactInterface* cinterface = m_contactInterfaces[ic].get();
        Eigen::SparseMatrix<real> sK(m_numDofs,m_numDofs);
        cinterface->compute_stiffness(sK);

        for (int k=0; k<sK.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(sK,k); it; ++it)
            {
                int rid = it.row();
                int cid = it.col();

                m_stiffness_matrix.coeffRef(rid,cid) += it.value();
            }
    }

    Eigen::VectorXd uload(m_numDofs);
    uload.setZero();

    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        zxNode* node =m_nodes[i].get();
        for(size_t j = 0; j < 3; j++)
        {
            int id = node->m_dof_id[j];
            if(id < -1)
            {
                id = -id - 2;
                uload[id] = node->rl[j] - node->rt[j];
            }
        }
    }

    m_Fd = m_stiffness_matrix * uload;

    for(size_t i = 0; i < m_numDofs; i++)
        if(m_is_prescribed_dofs[i])
            m_Fd[i] = 0.0;

    for (int k=0; k<m_stiffness_matrix.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(m_stiffness_matrix,k); it; ++it)
        {
            int rid = it.row();
            int cid = it.col();

            if(m_is_prescribed_dofs[rid] || m_is_prescribed_dofs[cid])
            {
                m_stiffness_matrix.coeffRef(rid,cid) = rid == cid? 1.0 : 0.0;
            }
        }

}

real zxALMSimulator::do_line_search(real s)
{

//    update_deformation(m_ui);
//    update_stress();
//    update_contact();

//    compute_residual(m_R1);

    //return 1.0;

    real smin = s;
    real a,A,B,D;
    real r0,r1,r;

    r0 = m_ui.dot(m_R0);

    real rmin = std::abs(r0);
    Eigen::VectorXd ul = m_ui;

    int n = 0;
    do
    {
        ul = m_ui;
        ul *= s;

        update_deformation(ul);
        update_stress();
        update_contact();

        compute_residual(m_R1);

        r1 = m_ui.dot(m_R1);

        if((n==0) || (std::abs(r1) < rmin))
        {
            smin = s;
            rmin = std::abs(r1);
        }

        if(std::abs(r1) < 1e-20) r = 0;
        else r = fabs(r1 / r0);

        if(r > 0.9)
        {
            a = r0/r1;
            A = 1 + a * (s-1);
            B = a * s*s;
            D = B * B - 4 * A * B;

            if(A == 0)
            {
                s = 1.0;
                break;
            }

            if(D >= 0)
            {
                s = (B + std::sqrt(D)) / (2 * A);
                if(s < 0) s = 0;
            }
            else
                s = 0.5 * B / A;
            ++n;
        }



    }while(n < 10 && r > 0.9);


    return s;

}


bool zxALMSimulator::alm_augment()
{
    bool converged = true;
    for(size_t ic = 0; ic < m_contactInterfaces.size(); ic++)
    {
        zxALMContactInterface* cinterface = m_contactInterfaces[ic].get();
        converged &= cinterface->alm_augment();
    }

    return converged;
}
