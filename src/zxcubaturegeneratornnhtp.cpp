#include "zxcubaturegeneratornnhtp.h"
#include "zxnnls_solver.h"

zxCubatureGeneratorNNHTP::zxCubatureGeneratorNNHTP(zxCubatureModel::Ptr cmodel)
{
    m_cubature_model = cmodel;

}

void zxCubatureGeneratorNNHTP::generateSamples(int numSample,double mag)
{
    m_cubature_model->generateSample(numSample,mag);
}

void zxCubatureGeneratorNNHTP::generateCubatures(int numC)
{
    int forceSize = m_cubature_model->get_Y_dim();
    int totalSamples = m_cubature_model->get_num_samples();
    int totalRows = forceSize * totalSamples;


    Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(m_cubature_model->get_sample_Y().data(),totalRows);
    int totalCandidates = m_cubature_model->get_num_elements();
    int maxKeyPoints = numC;
    double errorTolerance = 1e-3;

    Eigen::VectorXd gradient(totalCandidates),gradientNoW(totalCandidates);
    Eigen::VectorXd gradientS(totalCandidates),w(totalCandidates),oldW(totalCandidates);
    gradient.setZero();gradientNoW.setZero(); gradientS.setZero();
    w.setZero();oldW.setZero();

    std::vector<int> wIndex;

    int maxIteration = 100;

    zxNNLS_SOLVER nnls(totalRows,maxKeyPoints);
    nnls.maxIter() = 1000;
    Eigen::VectorXd bNNLS(totalRows);
    Eigen::VectorXd weightsNNLS(maxKeyPoints);
    weightsNNLS.setZero();
    Eigen::VectorXd A(totalRows * maxKeyPoints);
    double    bNorm = b.norm();
    double    rNorm = bNorm;
    double    relativeError = 1.0;
    wIndex = random_pick_candidates(wIndex,maxKeyPoints,totalCandidates,false);
    for(int i = 0; i < wIndex.size(); i++)
        w[wIndex[i]] = 1.0;

    std::map<int,Eigen::VectorXd> sampleColumns;

    for(int iter = 0; iter < maxIteration; iter++)
    {
        oldW = w;
        std::vector<int> candidateSubset = random_pick_candidates(wIndex,maxKeyPoints * 10,totalCandidates,true);


        for(size_t i = 0; i < candidateSubset.size(); i++)
        {
            if(sampleColumns.find(candidateSubset[i]) == sampleColumns.end())
                sampleColumns[candidateSubset[i]] = compute_element_Quantity(candidateSubset[i]);
        }

        Eigen::VectorXd Aw(totalRows);
        Aw.setZero();
        for(int i = 0; i < wIndex.size(); i++)
            Aw += w[wIndex[i]] * sampleColumns[wIndex[i]];
        Aw = b - Aw;

        gradient.setZero();
        for(unsigned int x = 0; x < candidateSubset.size(); x++){
            gradient[candidateSubset[x]] = 2 * (sampleColumns[candidateSubset[x]].dot(Aw));
        }

        gradientS = gradient;

        Eigen::VectorXd Aw_s(totalRows,1);
        Aw_s.setZero();
        for(unsigned int x = 0; x < candidateSubset.size(); x++){
            if(gradientS[candidateSubset[x]] > 0 || w[candidateSubset[x]] > 0)
                Aw_s += (gradientS[candidateSubset[x]] * sampleColumns[candidateSubset[x]]);
            else
                gradientS[candidateSubset[x]] = 0;
        }

        double stepSize = gradientS.squaredNorm() / Aw_s.squaredNorm();
        if(Aw_s.squaredNorm() == 0)
            continue;

        w += (stepSize * gradient);

        std::vector<int> newCandidates = project(w, maxKeyPoints);

        for(unsigned int x = 0; x < newCandidates.size(); x++){
            if(sampleColumns.find(newCandidates[x]) == sampleColumns.end()){

                sampleColumns[newCandidates[x]] = compute_element_Quantity(newCandidates[x]);
            }
            Eigen::VectorXd& column = sampleColumns[newCandidates[x]];
            for(int y = 0; y < totalRows; y++){
                A[x * totalRows + y] = column[y];
            }
        }
        for (int x = 0; x < totalRows; x++)
            bNNLS[x] = b[x];

        bool converged = nnls.solve(A.data(), newCandidates.size(), bNNLS.data(), weightsNNLS.data(), rNorm);

        relativeError = rNorm / bNorm;

        cout << "iteration " << iter << " relative error " << relativeError << endl;

        w.setZero();
        wIndex.clear();
        for(size_t x = 0; x < newCandidates.size(); x++){
            if(weightsNNLS[x] <= 1e-12)
                continue;

            int idx = newCandidates[x];
            w[idx] = weightsNNLS[x];
            wIndex.push_back(idx);
        }
        if(relativeError <= errorTolerance)
            break;

        if((w - oldW).squaredNorm() < 1e-6){
            cout << "picked same samples in two consecutive steps!!! stop" << endl;
            break;
        }

    }


    mKeyElements = wIndex;
    for(size_t i = 0; i < wIndex.size(); i++)
        mKeyWeights.push_back(w[wIndex[i]]);

}

std::vector<int> zxCubatureGeneratorNNHTP::random_pick_candidates(std::vector<int> &excludes, int candidatesPerTry, int totalCandidates, bool addExcludes)
{
    std::vector<bool> alreadyused(totalCandidates, false);
    for(unsigned int x = 0; x < excludes.size(); x++)
        alreadyused[excludes[x]] = true;

    MERSENNETWISTER& twister = m_twister;

    std::vector<int> candidates;

    if(addExcludes)
        candidates = excludes;

    candidatesPerTry = candidatesPerTry < (totalCandidates - candidates.size()) ? candidatesPerTry : (totalCandidates - candidates.size());

    while(candidates.size() < candidatesPerTry){
        int index = twister.randInt(totalCandidates - 1);
        if(!alreadyused[index]){
            candidates.push_back(index);
            alreadyused[index] = true;
        }
    }
    return candidates;
}

std::vector<int> zxCubatureGeneratorNNHTP::project(Eigen::VectorXd &w, int toKeep)
{
    struct orderByValueSmallerThan{
        bool operator()(std::pair<int, double> const& a, std::pair<int, double> const& b) const{
            return a.second < b.second;
        }
    };

    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, orderByValueSmallerThan> candidateHeap;
    for(int x = 0; x < w.size(); x++){
        if(w[x] > 0)
            candidateHeap.push(std::make_pair(x, w[x]));
    }
    w.setZero();
    std::vector<int> indices;
    toKeep = candidateHeap.size() < toKeep ? candidateHeap.size() : toKeep;
    for(int x = 0; x < toKeep; x++){
        std::pair<int, double> candidate = candidateHeap.top();

        indices.push_back(candidate.first);
        w[candidate.first] = candidate.second;
        candidateHeap.pop();
    }
    return indices;

}

Eigen::VectorXd zxCubatureGeneratorNNHTP::compute_element_Quantity(int id)
{
    return m_cubature_model->get_element_column(id);

}
