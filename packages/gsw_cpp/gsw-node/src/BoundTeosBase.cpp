//
// Created by blake on 9/26/24.
//

#include "BoundTeosBase.h"

Napi::Object BoundTeosBase::Init(Napi::Env env, Napi::Object exports)
{
    const auto func = DefineClass(
        env,
        "TeosBase",
        {
            InstanceMethod<&BoundTeosBase::gsw_z_from_p>(
                "gsw_z_from_p",
                static_cast<napi_property_attributes>(napi_writable | napi_configurable)
            ),
            InstanceMethod<&BoundTeosBase::gsw_sp_from_c>(
                "gsw_sp_from_c",
                static_cast<napi_property_attributes>(napi_writable | napi_configurable)
            ),
            InstanceMethod<&BoundTeosBase::gsw_sa_from_sp>(
                "gsw_sa_from_sp",
                static_cast<napi_property_attributes>(napi_writable | napi_configurable)
            ),
            InstanceMethod<&BoundTeosBase::gsw_ct_from_t>(
                "gsw_ct_from_t",
                static_cast<napi_property_attributes>(napi_writable | napi_configurable)
            ),
            InstanceMethod<&BoundTeosBase::gsw_depth_from_z>(
                "gsw_depth_from_z",
                static_cast<napi_property_attributes>(napi_writable | napi_configurable)
            )
        }
    );
    auto* constructor = new Napi::FunctionReference();
    *constructor = Persistent(func);
    (void)exports.Set("TeosBase", func);
    env.SetInstanceData<Napi::FunctionReference>(constructor);
    return exports;
}

BoundTeosBase::BoundTeosBase(const Napi::CallbackInfo& info): Napi::ObjectWrap<BoundTeosBase>(info)
{
}

Napi::Value BoundTeosBase::gsw_z_from_p(const Napi::CallbackInfo& info)
{
    const auto& env = info.Env();
    const auto& result = m_base.gsw_z_from_p(
        info[0].As<Napi::Number>().DoubleValue(),
        info[1].As<Napi::Number>().DoubleValue(),
        info[2].As<Napi::Number>().DoubleValue(),
        info[3].As<Napi::Number>().DoubleValue()
    );
    return Napi::Number::New(env, result);
}

Napi::Value BoundTeosBase::gsw_sp_from_c(const Napi::CallbackInfo& info)
{
    const auto& env = info.Env();
    const auto& result = m_base.gsw_sp_from_c(
        info[0].As<Napi::Number>().DoubleValue(),
        info[1].As<Napi::Number>().DoubleValue(),
        info[2].As<Napi::Number>().DoubleValue()
    );
    return Napi::Number::New(env, result);
}

Napi::Value BoundTeosBase::gsw_sa_from_sp(const Napi::CallbackInfo& info)
{
    const auto& env = info.Env();
    const auto& result = m_base.gsw_sa_from_sp(
        info[0].As<Napi::Number>().DoubleValue(),
        info[1].As<Napi::Number>().DoubleValue(),
        info[2].As<Napi::Number>().DoubleValue(),
        info[3].As<Napi::Number>().DoubleValue()
    );
    return Napi::Number::New(env, result);
}

Napi::Value BoundTeosBase::gsw_ct_from_t(const Napi::CallbackInfo& info)
{
    const auto& env = info.Env();
    const auto& result = m_base.gsw_ct_from_t(
        info[0].As<Napi::Number>().DoubleValue(),
        info[1].As<Napi::Number>().DoubleValue(),
        info[2].As<Napi::Number>().DoubleValue()
    );
    return Napi::Number::New(env, result);
}

Napi::Value BoundTeosBase::gsw_depth_from_z(const Napi::CallbackInfo& info)
{
    const auto& env = info.Env();
    const auto& result = m_base.gsw_depth_from_z(
        info[0].As<Napi::Number>().DoubleValue()
    );
    return Napi::Number::New(env, result);
}
