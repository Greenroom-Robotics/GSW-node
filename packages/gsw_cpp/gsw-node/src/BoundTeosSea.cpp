//
// Created by blake on 9/26/24.
//

#include "BoundTeosSea.h"

Napi::Object BoundTeosSea::Init(Napi::Env env, Napi::Object exports)
{
    const auto func = DefineClass(
        env,
        "TeosSea",
        {
            InstanceMethod<&BoundTeosSea::gsw_c_from_sp>(
                "gsw_c_from_sp",
                static_cast<napi_property_attributes>(napi_writable | napi_configurable)
            ),
            InstanceMethod<&BoundTeosSea::gsw_sound_speed>(
                "gsw_sound_speed",
                static_cast<napi_property_attributes>(napi_writable | napi_configurable)
            ),
        }
    );
    auto* constructor = new Napi::FunctionReference();
    *constructor = Persistent(func);
    (void)exports.Set("TeosSea", func);
    env.SetInstanceData<Napi::FunctionReference>(constructor);
    return exports;
}

BoundTeosSea::BoundTeosSea(const Napi::CallbackInfo& info): Napi::ObjectWrap<BoundTeosSea>(info)
{
}

Napi::Value BoundTeosSea::gsw_c_from_sp(const Napi::CallbackInfo& info)
{
    const auto& env = info.Env();
    const auto& result = m_sea.gsw_c_from_sp(
        info[0].As<Napi::Number>().DoubleValue(),
        info[1].As<Napi::Number>().DoubleValue(),
        info[2].As<Napi::Number>().DoubleValue()
    );
    return Napi::Number::New(env, result);
}

Napi::Value BoundTeosSea::gsw_sound_speed(const Napi::CallbackInfo& info)
{
    const auto& env = info.Env();
    const auto& result = m_sea.gsw_sound_speed(
        info[0].As<Napi::Number>().DoubleValue(),
        info[1].As<Napi::Number>().DoubleValue(),
        info[2].As<Napi::Number>().DoubleValue()
    );
    return Napi::Number::New(env, result);
}
