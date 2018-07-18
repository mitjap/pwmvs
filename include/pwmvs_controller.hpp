#pragma once

#include "types.hpp"
#include "workspace.hpp"

class AbstractProgress;

template <class T>
class Controller
{
public:
    Controller(const std::shared_ptr<Workspace> &workspace, const typename T::Options &options = T::Options());
    bool run(bool geometric, AbstractProgress *progress = nullptr);

private:
    bool runInternal(bool geometric, AbstractProgress *progress);

public:
    std::shared_ptr<Workspace> workspace;
    typename T::Options pwmvs_options;
};
