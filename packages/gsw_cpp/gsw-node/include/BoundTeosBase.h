//
// Created by blake on 9/26/24.
//

#ifndef BOUNDTEOSBASE_H
#define BOUNDTEOSBASE_H

#include <napi.h>

#include <TeosBase.h>

class BoundTeosBase : public Napi::ObjectWrap<BoundTeosBase>
{
public:
    static Napi::Object Init(Napi::Env env, Napi::Object exports);

    explicit BoundTeosBase(const Napi::CallbackInfo& info);

    Napi::Value gsw_z_from_p(const Napi::CallbackInfo& info);

    Napi::Value gsw_sp_from_c(const Napi::CallbackInfo& info);

    Napi::Value gsw_sa_from_sp(const Napi::CallbackInfo& info);

    Napi::Value gsw_ct_from_t(const Napi::CallbackInfo& info);

    Napi::Value gsw_depth_from_z(const Napi::CallbackInfo& info);

private:
    TeosBase m_base;
};


#endif //BOUNDTEOSBASE_H
