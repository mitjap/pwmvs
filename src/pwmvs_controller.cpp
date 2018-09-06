#include "pwmvs_controller.hpp"
#include "pwmvs.hpp"

#include "view_utils.hpp"
#include "progress.hpp"

#include <stlplus3/filesystemSimplified/file_system.hpp>


template <class T>
Controller<T>::Controller(const std::shared_ptr<Workspace> &workspace, const typename T::Options &pwmvs_options)
    : workspace(workspace)
    , pwmvs_options(pwmvs_options)
{
    // nothing to do
}

template <class T>
bool Controller<T>::run(bool geometric, AbstractProgress *progress)
{
    if (!stlplus::folder_exists(workspace->work_path) && !stlplus::folder_create(workspace->work_path))
    {
        std::cerr << "Can not create working folder." << std::endl;
        return false;
    }

    if (progress) {
        progress->configure((geometric ? 2 : 1) * workspace->view_data.size());
    }

    if (!runInternal(false, progress))
        return false;
    if (geometric && !runInternal(true, progress))
        return false;

    return true;
}

template <class T>
bool Controller<T>::runInternal(bool geometric, AbstractProgress *progress)
{
    for (int ref_id = 0; ref_id < workspace->view_data.size(); ++ref_id)
    {
        ProgressIncrementor inc(progress);

        const ViewData &ref_view = workspace->view_data[ref_id];
        if (isComplete(*workspace, ref_id, geometric))
            continue;

        FloatT min_depth, max_depth;
        if (!workspace->getMinMaxDepth(ref_id, min_depth, max_depth))
            continue;

        std::vector<int> src_view_ids;
        if (!workspace->getSrcViewIds(ref_id, src_view_ids))
            continue;

        std::cout << "Source views:";
        for (int src_view_id : src_view_ids)
            std::cout << " " << src_view_id;
        std::cout << std::endl;

        std::shared_ptr<RefView> ref = createRefView(ref_view);
        if (!loadFromFiles(*workspace, ref_id, false, true, false, geometric, geometric, *ref)) {
            continue;
        }

        std::vector<std::shared_ptr<SrcView>> srcs;
        for (int src_id : src_view_ids)
        {
            if (srcs.size() >= pwmvs_options.max_sources && pwmvs_options.max_sources > 0)
                break;

            const ViewData &src_view = workspace->view_data[src_id];

            std::shared_ptr<SrcView> src = createSrcView(ref_view, src_view);
            if (!loadFromFiles(*workspace, src_id, false, true, false, geometric, geometric, *src)) {
                continue;
            }

            srcs.push_back(src);
        }

        if (!srcs.empty())
        {
            typename T::Options options = pwmvs_options;
            options.min_depth = min_depth;
            options.max_depth = max_depth;
            options.geometric_consistency_term &= geometric;
            options.filter_geometric_consistency &= geometric;
            options.filter_min_num_consistent = std::min(options.filter_min_num_consistent, static_cast<decltype(options.filter_min_num_consistent)>(srcs.size()));
            T pwmvs(ref, srcs, options);
            if (!pwmvs.run())
                return false;

            saveToFiles(*workspace, ref_id, true, true, geometric, *ref);
            //debugToFiles(workspace, ref_id, true, true, geometric, *ref, min_depth, max_depth);
        }
    }

    return true;
}

template class Controller<PWMVS>;
