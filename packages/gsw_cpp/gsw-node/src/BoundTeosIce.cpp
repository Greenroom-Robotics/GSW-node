//
// Created by blake on 9/26/24.
//

#include "BoundTeosIce.h"

Napi::Object BoundTeosIce::Init(Napi::Env env, Napi::Object exports)
{
    const auto func = DefineClass(
        env,
        "TeosIce",
        {
            InstanceMethod<&BoundTeosIce::gsw_cp_ice>(
                "gsw_cp_ice",
                static_cast<napi_property_attributes>(napi_writable | napi_configurable)
            ),
        }
    );
    auto* constructor = new Napi::FunctionReference();
    *constructor = Persistent(func);
    (void)exports.Set("TeosIce", func);
    env.SetInstanceData<Napi::FunctionReference>(constructor);
    return exports;
}

BoundTeosIce::BoundTeosIce(const Napi::CallbackInfo& info): Napi::ObjectWrap<BoundTeosIce>(info)
{
}

Napi::Value BoundTeosIce::gsw_cp_ice(const Napi::CallbackInfo& info)
{
    const auto& env = info.Env();
    const auto& result = m_ice.gsw_cp_ice(
        info[0].As<Napi::Number>().DoubleValue(),
        info[1].As<Napi::Number>().DoubleValue()
    );
    return Napi::Number::New(env, result);
}
