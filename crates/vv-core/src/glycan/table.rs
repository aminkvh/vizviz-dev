//! Residue-name -> monosaccharide table, transcribed from the
//! `3D-SNFG.tcl` script (version 1, D.F. Thieker & J.A. Hadden, 2016): every
//! `_common`/`_charmm`/`_glycam` name list in its `SNFG` namespace, paired
//! with the shape/color its `snfg-detect` assigns that group. Entry order
//! matches the source's `elseif` chain, which matters for the one name
//! (`4YS`) two groups share: GlcNAc (tested first) wins there, as it does
//! in the source. See `super` for the citation.
use super::{MonoEntry, Shape, SnfgColor};

pub(super) const MONO_TABLE: &[MonoEntry] = &[
    MonoEntry {
        key: "Glc",
        label: "Glucose (blue sphere)",
        shape: Shape::Sphere,
        color1: SnfgColor::Blue,
        color2: SnfgColor::Blue,
        names: &[
            "GLC", "MAL", "BGC", "AGLC", "BGLC", "0GA", "0GB", "1GA", "1GB", "2GA", "2GB", "3GA",
            "3GB", "4GA", "4GB", "6GA", "6GB", "ZGA", "ZGB", "YGA", "YGB", "XGA", "XGB", "WGA",
            "WGB", "VGA", "VGB", "UGA", "UGB", "TGA", "TGB", "SGA", "SGB", "RGA", "RGB", "QGA",
            "QGB", "PGA", "PGB", "0gA", "0gB", "1gA", "1gB", "2gA", "2gB", "3gA", "3gB", "4gA",
            "4gB", "6gA", "6gB", "ZgA", "ZgB", "YgA", "YgB", "XgA", "XgB", "WgA", "WgB", "VgA",
            "VgB", "UgA", "UgB", "TgA", "TgB", "SgA", "SgB", "RgA", "RgB", "QgA", "QgB", "PgA",
            "PgB",
        ],
    },
    MonoEntry {
        key: "Man",
        label: "Mannose (green sphere)",
        shape: Shape::Sphere,
        color1: SnfgColor::Green,
        color2: SnfgColor::Green,
        names: &[
            "MAN", "BMA", "AMAN", "BMAN", "0MA", "0MB", "1MA", "1MB", "2MA", "2MB", "3MA", "3MB",
            "4MA", "4MB", "6MA", "6MB", "ZMA", "ZMB", "YMA", "YMB", "XMA", "XMB", "WMA", "WMB",
            "VMA", "VMB", "UMA", "UMB", "TMA", "TMB", "SMA", "SMB", "RMA", "RMB", "QMA", "QMB",
            "PMA", "PMB", "0mA", "0mB", "1mA", "1mB", "2mA", "2mB", "3mA", "3mB", "4mA", "4mB",
            "6mA", "6mB", "ZmA", "ZmB", "YmA", "YmB", "XmA", "XmB", "WmA", "WmB", "VmA", "VmB",
            "UmA", "UmB", "TmA", "TmB", "SmA", "SmB", "RmA", "RmB", "QmA", "QmB", "PmA", "PmB",
        ],
    },
    MonoEntry {
        key: "Gal",
        label: "Galactose (yellow sphere)",
        shape: Shape::Sphere,
        color1: SnfgColor::Yellow,
        color2: SnfgColor::Yellow,
        names: &[
            "GAL", "GLA", "AGAL", "BGAL", "0LA", "0LB", "1LA", "1LB", "2LA", "2LB", "3LA", "3LB",
            "4LA", "4LB", "6LA", "6LB", "ZLA", "ZLB", "YLA", "YLB", "XLA", "XLB", "WLA", "WLB",
            "VLA", "VLB", "ULA", "ULB", "TLA", "TLB", "SLA", "SLB", "RLA", "RLB", "QLA", "QLB",
            "PLA", "PLB", "0lA", "0lB", "1lA", "1lB", "2lA", "2lB", "3lA", "3lB", "4lA", "4lB",
            "6lA", "6lB", "ZlA", "ZlB", "YlA", "YlB", "XlA", "XlB", "WlA", "WlB", "VlA", "VlB",
            "UlA", "UlB", "TlA", "TlB", "SlA", "SlB", "RlA", "RlB", "QlA", "QlB", "PlA", "PlB",
        ],
    },
    MonoEntry {
        key: "Gul",
        label: "Gulose (orange sphere)",
        shape: Shape::Sphere,
        color1: SnfgColor::Orange,
        color2: SnfgColor::Orange,
        names: &[
            "GUL", "GUP", "GL0", "AGUL", "BGUL", "0KA", "0KB", "1KA", "1KB", "2KA", "2KB", "3KA",
            "3KB", "4KA", "4KB", "6KA", "6KB", "ZKA", "ZKB", "YKA", "YKB", "XKA", "XKB", "WKA",
            "WKB", "VKA", "VKB", "UKA", "UKB", "TKA", "TKB", "SKA", "SKB", "RKA", "RKB", "QKA",
            "QKB", "PKA", "PKB", "0kA", "0kB", "1kA", "1kB", "2kA", "2kB", "3kA", "3kB", "4kA",
            "4kB", "6kA", "6kB", "ZkA", "ZkB", "YkA", "YkB", "XkA", "XkB", "WkA", "WkB", "VkA",
            "VkB", "UkA", "UkB", "TkA", "TkB", "SkA", "SkB", "RkA", "RkB", "QkA", "QkB", "PkA",
            "PkB",
        ],
    },
    MonoEntry {
        key: "Alt",
        label: "Altrose (pink sphere)",
        shape: Shape::Sphere,
        color1: SnfgColor::Pink,
        color2: SnfgColor::Pink,
        names: &[
            "ALT", "AALT", "BALT", "0EA", "0EB", "1EA", "1EB", "2EA", "2EB", "3EA", "3EB", "4EA",
            "4EB", "6EA", "6EB", "ZEA", "ZEB", "YEA", "YEB", "XEA", "XEB", "WEA", "WEB", "VEA",
            "VEB", "UEA", "UEB", "TEA", "TEB", "SEA", "SEB", "REA", "REB", "QEA", "QEB", "PEA",
            "PEB", "0eA", "0eB", "1eA", "1eB", "2eA", "2eB", "3eA", "3eB", "4eA", "4eB", "6eA",
            "6eB", "ZeA", "ZeB", "YeA", "YeB", "XeA", "XeB", "WeA", "WeB", "VeA", "VeB", "UeA",
            "UeB", "TeA", "TeB", "SeA", "SeB", "ReA", "ReB", "QeA", "QeB", "PeA", "PeB",
        ],
    },
    MonoEntry {
        key: "All",
        label: "Allose (purple sphere)",
        shape: Shape::Sphere,
        color1: SnfgColor::Purple,
        color2: SnfgColor::Purple,
        names: &[
            "ALL", "WOO", "AALL", "BALL", "0NA", "0NB", "1NA", "1NB", "2NA", "2NB", "3NA", "3NB",
            "4NA", "4NB", "6NA", "6NB", "ZNA", "ZNB", "YNA", "YNB", "XNA", "XNB", "WNA", "WNB",
            "VNA", "VNB", "UNA", "UNB", "TNA", "TNB", "SNA", "SNB", "RNA", "RNB", "QNA", "QNB",
            "PNA", "PNB", "0nA", "0nB", "1nA", "1nB", "2nA", "2nB", "3nA", "3nB", "4nA", "4nB",
            "6nA", "6nB", "ZnA", "ZnB", "YnA", "YnB", "XnA", "XnB", "WnA", "WnB", "VnA", "VnB",
            "UnA", "UnB", "TnA", "TnB", "SnA", "SnB", "RnA", "RnB", "QnA", "QnB", "PnA", "PnB",
        ],
    },
    MonoEntry {
        key: "Tal",
        label: "Talose (light blue sphere)",
        shape: Shape::Sphere,
        color1: SnfgColor::LightBlue,
        color2: SnfgColor::LightBlue,
        names: &[
            "TAL", "ATAL", "BTAL", "0TA", "0TB", "1TA", "1TB", "2TA", "2TB", "3TA", "3TB", "4TA",
            "4TB", "6TA", "6TB", "ZTA", "ZTB", "YTA", "YTB", "XTA", "XTB", "WTA", "WTB", "VTA",
            "VTB", "UTA", "UTB", "TTA", "TTB", "STA", "STB", "RTA", "RTB", "QTA", "QTB", "PTA",
            "PTB", "0tA", "0tB", "1tA", "1tB", "2tA", "2tB", "3tA", "3tB", "4tA", "4tB", "6tA",
            "6tB", "ZtA", "ZtB", "YtA", "YtB", "XtA", "XtB", "WtA", "WtB", "VtA", "VtB", "UtA",
            "UtB", "TtA", "TtB", "StA", "StB", "RtA", "RtB", "QtA", "QtB", "PtA", "PtB",
        ],
    },
    MonoEntry {
        key: "Ido",
        label: "Idose (brown sphere)",
        shape: Shape::Sphere,
        color1: SnfgColor::Brown,
        color2: SnfgColor::Brown,
        names: &["IDO", "AIDO", "BIDO"],
    },
    MonoEntry {
        key: "GlcNAc",
        label: "N-acetyl-glucosamine (blue cube)",
        shape: Shape::Cube,
        color1: SnfgColor::Blue,
        color2: SnfgColor::Blue,
        names: &[
            "NAG", "4YS", "SGN", "BGLN", "NDG", "AGLCNA", "BGLCNA", "BGLCN0", "0YA", "0YB", "1YA",
            "1YB", "3YA", "3YB", "4YA", "4YB", "6YA", "6YB", "WYA", "WYB", "VYA", "VYB", "UYA",
            "UYB", "QYA", "QYB", "0yA", "0yB", "1yA", "1yB", "3yA", "3yB", "4yA", "4yB", "6yA",
            "6yB", "WyA", "WyB", "VyA", "VyB", "UyA", "UyB", "QyA", "QyB", "XYY", "UYY", "VYY",
        ],
    },
    MonoEntry {
        key: "ManNAc",
        label: "N-acetyl-mannosamine (green cube)",
        shape: Shape::Cube,
        color1: SnfgColor::Green,
        color2: SnfgColor::Green,
        names: &[
            "0WA", "0WB", "1WA", "1WB", "3WA", "3WB", "4WA", "4WB", "6WA", "6WB", "WWA", "WWB",
            "VWA", "VWB", "UWA", "UWB", "QWA", "QWB", "0wA", "0wB", "1wA", "1wB", "3wA", "3wB",
            "4wA", "4wB", "6wA", "6wB", "WwA", "WwB", "VwA", "VwB", "UwA", "UwB", "QwA", "QwB",
        ],
    },
    MonoEntry {
        key: "GalNAc",
        label: "N-acetyl-galactosamine (yellow cube)",
        shape: Shape::Cube,
        color1: SnfgColor::Yellow,
        color2: SnfgColor::Yellow,
        names: &[
            "NGA", "AGALNA", "BGALNA", "0VA", "0VB", "1VA", "1VB", "3VA", "3VB", "4VA", "4VB",
            "6VA", "6VB", "WVA", "WVB", "VVA", "VVB", "UVA", "UVB", "QVA", "QVB", "0vA", "0vB",
            "1vA", "1vB", "3vA", "3vB", "4vA", "4vB", "6vA", "6vB", "WvA", "WvB", "VvA", "VvB",
            "UvA", "UvB", "QvA", "QvB",
            "A2G", // CCD: 2-acetamido-2-deoxy-alpha-D-galactopyranose
        ],
    },
    MonoEntry {
        key: "GulNAc",
        label: "N-Acetyl-D-gulosamine (orange cube)",
        shape: Shape::Cube,
        color1: SnfgColor::Orange,
        color2: SnfgColor::Orange,
        names: &[],
    },
    MonoEntry {
        key: "AltNAc",
        label: "N-Acetyl-L-altrosamine (pink cube)",
        shape: Shape::Cube,
        color1: SnfgColor::Pink,
        color2: SnfgColor::Pink,
        names: &[],
    },
    MonoEntry {
        key: "AllNAc",
        label: "N-Acetyl-D-allosamine (purple cube)",
        shape: Shape::Cube,
        color1: SnfgColor::Purple,
        color2: SnfgColor::Purple,
        names: &[],
    },
    MonoEntry {
        key: "TalNAc",
        label: "N-Acetyl-D-talosamine (light blue cube)",
        shape: Shape::Cube,
        color1: SnfgColor::LightBlue,
        color2: SnfgColor::LightBlue,
        names: &[],
    },
    MonoEntry {
        key: "IdoNAc",
        label: "N-Acetyl-L-idosamine (brown cube)",
        shape: Shape::Cube,
        color1: SnfgColor::Brown,
        color2: SnfgColor::Brown,
        names: &[],
    },
    MonoEntry {
        key: "GlcN",
        label: "Glucosamine (white\\blue cube)",
        shape: Shape::Cube,
        color1: SnfgColor::White,
        color2: SnfgColor::Blue,
        names: &[
            "GCS", "0YN", "0Yn", "0YNP", "0YnP", "0YS", "0Ys", "3YS", "3Ys", "4YS", "4Ys", "6YS",
            "6Ys", "QYS", "QYs", "UYS", "UYs", "VYS", "VYs", "WYS", "WYs", "0yS", "0ys", "3yS",
            "3ys", "4yS", "4ys", "PA1", // CCD: 2-amino-2-deoxy-alpha-D-glucopyranose
        ],
    },
    MonoEntry {
        key: "ManN",
        label: "D-Mannosamine (white\\green cube)",
        shape: Shape::Cube,
        color1: SnfgColor::White,
        color2: SnfgColor::Green,
        names: &[],
    },
    MonoEntry {
        key: "GalN",
        label: "D-Galactosamine (white\\yellow cube)",
        shape: Shape::Cube,
        color1: SnfgColor::White,
        color2: SnfgColor::Yellow,
        names: &["X6X"], // CCD: 2-amino-2-deoxy-alpha-D-galactopyranose
    },
    MonoEntry {
        key: "GulN",
        label: "D-Gulosamine (white\\orange cube)",
        shape: Shape::Cube,
        color1: SnfgColor::White,
        color2: SnfgColor::Orange,
        names: &[],
    },
    MonoEntry {
        key: "AltN",
        label: "L-Altrosamine (white\\pink cube)",
        shape: Shape::Cube,
        color1: SnfgColor::White,
        color2: SnfgColor::Pink,
        names: &[],
    },
    MonoEntry {
        key: "AllN",
        label: "D-Allosamine (white\\purple cube)",
        shape: Shape::Cube,
        color1: SnfgColor::White,
        color2: SnfgColor::Purple,
        names: &[],
    },
    MonoEntry {
        key: "TalN",
        label: "D-Talosamine (white\\light blue cube)",
        shape: Shape::Cube,
        color1: SnfgColor::White,
        color2: SnfgColor::LightBlue,
        names: &[],
    },
    MonoEntry {
        key: "IdoN",
        label: "L-Idosamine (white\\brown cube)",
        shape: Shape::Cube,
        color1: SnfgColor::White,
        color2: SnfgColor::Brown,
        names: &[],
    },
    MonoEntry {
        key: "GlcA",
        label: "Glucuronic acid (blue-white diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::Blue,
        color2: SnfgColor::White,
        names: &[
            "GCU", "AGLCA", "BGLCA", "BGLCA0", "0ZA", "0ZB", "1ZA", "1ZB", "2ZA", "2ZB", "3ZA",
            "3ZB", "4ZA", "4ZB", "ZZA", "ZZB", "YZA", "YZB", "WZA", "WZB", "TZA", "TZB", "0zA",
            "0zB", "1zA", "1zB", "2zA", "2zB", "3zA", "3zB", "4zA", "4zB", "ZzA", "ZzB", "YzA",
            "YzB", "WzA", "WzB", "TzA", "TzB", "0ZBP",
            "BDP", // CCD: beta-D-glucopyranuronic acid
        ],
    },
    MonoEntry {
        key: "ManA",
        label: "Mannuronic acid (green-white diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::Green,
        color2: SnfgColor::White,
        names: &["MAV", "BEM"],
    },
    MonoEntry {
        key: "GalA",
        label: "Galacturonic acid (yellow-white diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::Yellow,
        color2: SnfgColor::White,
        names: &[
            "ADA", "0OA", "0OB", "1OA", "1OB", "2OA", "2OB", "3OA", "3OB", "4OA", "4OB", "ZOA",
            "ZOB", "YOA", "YOB", "WOA", "WOB", "TOA", "TOB", "0oA", "0oB", "1oA", "1oB", "2oA",
            "2oB", "3oA", "3oB", "4oA", "4oB", "ZoA", "ZoB", "YoA", "YoB", "WoA", "WoB", "ToA",
            "ToB", "GTR", // CCD: beta-D-galactopyranuronic acid
        ],
    },
    MonoEntry {
        key: "AltA",
        label: "L-Altruronic acid (white-pink diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::White,
        color2: SnfgColor::Pink,
        names: &[],
    },
    MonoEntry {
        key: "AllA",
        label: "D-Alluronic acid (purple-white diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::Purple,
        color2: SnfgColor::White,
        names: &[],
    },
    MonoEntry {
        key: "TalA",
        label: "D-Taluronic acid (light blue-white diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::LightBlue,
        color2: SnfgColor::White,
        names: &[],
    },
    MonoEntry {
        key: "GulA",
        label: "Guluronic acid (orange-white diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::Orange,
        color2: SnfgColor::White,
        names: &["LGU"],
    },
    MonoEntry {
        key: "IdoA",
        label: "Iduronic acid (white-brown diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::White,
        color2: SnfgColor::Brown,
        names: &[
            "IDS", "AIDOA", "BIDOA", "0UA", "0UB", "1UA", "1UB", "2UA", "2UB", "3UA", "3UB", "4UA",
            "4UB", "ZUA", "ZUB", "YUA", "YUB", "WUA", "WUB", "TUA", "TUB", "0uA", "0uB", "1uA",
            "1uB", "2uA", "2uB", "3uA", "3uB", "4uA", "4uB", "ZuA", "ZuB", "YuA", "YuB", "WuA",
            "WuB", "TuA", "TuB", "YuAP",
            "IDR", // CCD: alpha-L-idopyranuronic acid (heparin/HS)
        ],
    },
    MonoEntry {
        key: "Qui",
        label: "Quinovose (blue cone)",
        shape: Shape::Cone,
        color1: SnfgColor::Blue,
        color2: SnfgColor::Blue,
        names: &[
            "QUI", "0QA", "0QB", "1QA", "1QB", "2QA", "2QB", "3QA", "3QB", "4QA", "4QB", "ZQA",
            "ZQB", "YQA", "YQB", "WQA", "WQB", "TQA", "TQB", "0qA", "0qB", "1qA", "1qB", "2qA",
            "2qB", "3qA", "3qB", "4qA", "4qB", "ZqA", "ZqB", "YqA", "YqB", "WqA", "WqB", "TqA",
            "TqB",
        ],
    },
    MonoEntry {
        key: "Rha",
        label: "Rhamnose (green cone)",
        shape: Shape::Cone,
        color1: SnfgColor::Green,
        color2: SnfgColor::Green,
        names: &[
            "RAM", "ARHM", "BRHM", "0HA", "0HB", "1HA", "1HB", "2HA", "2HB", "3HA", "3HB", "4HA",
            "4HB", "ZHA", "ZHB", "YHA", "YHB", "WHA", "WHB", "THA", "THB", "0hA", "0hB", "1hA",
            "1hB", "2hA", "2hB", "3hA", "3hB", "4hA", "4hB", "ZhA", "ZhB", "YhA", "YhB", "WhA",
            "WhB", "ThA", "ThB", "RM4", // CCD: beta-L-rhamnopyranose
        ],
    },
    MonoEntry {
        key: "Fuc",
        label: "Fucose (red cone)",
        shape: Shape::Cone,
        color1: SnfgColor::Red,
        color2: SnfgColor::Red,
        names: &[
            "FUC", "FUL", "AFUC", "BFUC", "0FA", "0FB", "1FA", "1FB", "2FA", "2FB", "3FA", "3FB",
            "4FA", "4FB", "ZFA", "ZFB", "YFA", "YFB", "WFA", "WFB", "TFA", "TFB", "0fA", "0fB",
            "1fA", "1fB", "2fA", "2fB", "3fA", "3fB", "4fA", "4fB", "ZfA", "ZfB", "YfA", "YfB",
            "WfA", "WfB", "TfA", "TfB",
        ],
    },
    MonoEntry {
        key: "6dGul",
        label: "6-Deoxy-D-gulose (orange cone)",
        shape: Shape::Cone,
        color1: SnfgColor::Orange,
        color2: SnfgColor::Orange,
        names: &[],
    },
    MonoEntry {
        key: "6dAlt",
        label: "6-Deoxy-L-altrose (pink cone)",
        shape: Shape::Cone,
        color1: SnfgColor::Pink,
        color2: SnfgColor::Pink,
        names: &[],
    },
    MonoEntry {
        key: "6dTal",
        label: "6-Deoxy-D-talose (light blue cone)",
        shape: Shape::Cone,
        color1: SnfgColor::LightBlue,
        color2: SnfgColor::LightBlue,
        names: &[],
    },
    MonoEntry {
        key: "QuiNAc",
        label: "N-Acetyl-D-quinovosamine (white-blue cone)",
        shape: Shape::Cone,
        color1: SnfgColor::White,
        color2: SnfgColor::Blue,
        names: &[],
    },
    MonoEntry {
        key: "RhaNAc",
        label: "N-Acetyl-L-rhamnosamine (white-green cone)",
        shape: Shape::Cone,
        color1: SnfgColor::White,
        color2: SnfgColor::Green,
        names: &[],
    },
    MonoEntry {
        key: "6dAltNAc",
        label: "N-Acetyl-6-deoxy-L-altrosamine (white-pink cone)",
        shape: Shape::Cone,
        color1: SnfgColor::White,
        color2: SnfgColor::Pink,
        names: &[],
    },
    MonoEntry {
        key: "6dTalNAc",
        label: "N-Acetyl-6-deoxy-D-talosamine (white-light blue cone)",
        shape: Shape::Cone,
        color1: SnfgColor::White,
        color2: SnfgColor::LightBlue,
        names: &[],
    },
    MonoEntry {
        key: "FucNAc",
        label: "N-Acetyl-L-fucosamine (white-red cone)",
        shape: Shape::Cone,
        color1: SnfgColor::White,
        color2: SnfgColor::Red,
        names: &[],
    },
    MonoEntry {
        key: "Oli",
        label: "Olivose (blue rectangle)",
        shape: Shape::Rectangle,
        color1: SnfgColor::Blue,
        color2: SnfgColor::Blue,
        names: &["OLI"],
    },
    MonoEntry {
        key: "Tyv",
        label: "Tyvelose (green rectangle)",
        shape: Shape::Rectangle,
        color1: SnfgColor::Green,
        color2: SnfgColor::Green,
        names: &[
            "TYV", "0TV", "0Tv", "1TV", "1Tv", "2TV", "2Tv", "4TV", "4Tv", "YTV", "YTv", "0tV",
            "0tv", "1tV", "1tv", "2tV", "2tv", "4tV", "4tv", "YtV", "Ytv",
        ],
    },
    MonoEntry {
        key: "Abe",
        label: "Abequose (orange rectangle)",
        shape: Shape::Rectangle,
        color1: SnfgColor::Orange,
        color2: SnfgColor::Orange,
        names: &[
            "ABE", "0AE", "2AE", "4AE", "YGa", "0AF", "2AF", "4AF", "YAF",
        ],
    },
    MonoEntry {
        key: "Par",
        label: "Paratose (pink rectangle)",
        shape: Shape::Rectangle,
        color1: SnfgColor::Pink,
        color2: SnfgColor::Pink,
        names: &["PAR"],
    },
    MonoEntry {
        key: "Dig",
        label: "Digitoxose (purple rectangle)",
        shape: Shape::Rectangle,
        color1: SnfgColor::Purple,
        color2: SnfgColor::Purple,
        names: &["DIG"],
    },
    MonoEntry {
        key: "Col",
        label: "Colitose (light blue rectangle)",
        shape: Shape::Rectangle,
        color1: SnfgColor::LightBlue,
        color2: SnfgColor::LightBlue,
        names: &["COL"],
    },
    MonoEntry {
        key: "Ara",
        label: "Arabinose (green star)",
        shape: Shape::Star,
        color1: SnfgColor::Green,
        color2: SnfgColor::Green,
        names: &[
            "ARA", "AHR", "AARB", "BARB", "0AA", "0AB", "1AA", "1AB", "2AA", "2AB", "3AA", "3AB",
            "4AA", "4AB", "ZAA", "ZAB", "YAA", "YAB", "WAA", "WAB", "TAA", "TAB", "0AD", "0AU",
            "1AD", "1AU", "2AD", "2AU", "3AD", "3AU", "5AD", "5AU", "ZAD", "ZAU", "0aA", "0aB",
            "1aA", "1aB", "2aA", "2aB", "3aA", "3aB", "4aA", "4aB", "ZaA", "ZaB", "YaA", "YaB",
            "WaA", "WaB", "TaA", "TaB", "0aD", "0aU", "1aD", "1aU", "2aD", "2aU", "3aD", "3aU",
            "5aD", "5aU", "ZaD", "ZaU", "ARB", // CCD: beta-L-arabinopyranose
        ],
    },
    MonoEntry {
        key: "Lyx",
        label: "Lyxose (yellow star)",
        shape: Shape::Star,
        color1: SnfgColor::Yellow,
        color2: SnfgColor::Yellow,
        names: &[
            "LYX", "ALYF", "BLYF", "0DA", "0DB", "1DA", "1DB", "2DA", "2DB", "3DA", "3DB", "4DA",
            "4DB", "ZDA", "ZDB", "YDA", "YDB", "WDA", "WDB", "TDA", "TDB", "0DD", "0DU", "1DD",
            "1DU", "2DD", "2DU", "3DD", "3DU", "5DD", "5DU", "ZDD", "ZDU", "0dA", "0dB", "1dA",
            "1dB", "2dA", "2dB", "3dA", "3dB", "4dA", "4dB", "ZdA", "ZdB", "YdA", "YdB", "WdA",
            "WdB", "TdA", "TdB", "0dD", "0dU", "1dD", "1dU", "2dD", "2dU", "3dD", "3dU", "5dD",
            "5dU", "ZdD", "ZdU",
        ],
    },
    MonoEntry {
        key: "Xyl",
        label: "Xylose (orange star)",
        shape: Shape::Star,
        color1: SnfgColor::Orange,
        color2: SnfgColor::Orange,
        names: &[
            "XYL", "XYS", "LXC", "XYP", "AXYL", "BXYL", "AXYF", "BXYF", "0XA", "0XB", "1XA", "1XB",
            "2XA", "2XB", "3XA", "3XB", "4XA", "4XB", "ZXA", "ZXB", "YXA", "YXB", "WXA", "WXB",
            "TXA", "TXB", "0XD", "0XU", "1XD", "1XU", "2XD", "2XU", "3XD", "3XU", "5XD", "5XU",
            "ZXD", "ZXU", "0xA", "0xB", "1xA", "1xB", "2xA", "2xB", "3xA", "3xB", "4xA", "4xB",
            "ZxA", "ZxB", "YxA", "YxB", "WxA", "WxB", "TxA", "TxB", "0xD", "0xU", "1xD", "1xU",
            "2xD", "2xU", "3xD", "3xU", "5xD", "5xU", "ZxD", "ZxU",
        ],
    },
    MonoEntry {
        key: "Rib",
        label: "Ribose (pink star)",
        shape: Shape::Star,
        color1: SnfgColor::Pink,
        color2: SnfgColor::Pink,
        names: &[
            "RIB", "ARIB", "BRIB", "0RA", "0RB", "1RA", "1RB", "2RA", "2RB", "3RA", "3RB", "4RA",
            "4RB", "ZRA", "ZRB", "YRA", "YRB", "WRA", "WRB", "TRA", "TRB", "0RD", "0RU", "1RD",
            "1RU", "2RD", "2RU", "3RD", "3RU", "5RD", "5RU", "ZRD", "ZRU", "0rA", "0rB", "1rA",
            "1rB", "2rA", "2rB", "3rA", "3rB", "4rA", "4rB", "ZrA", "ZrB", "YrA", "YrB", "WrA",
            "WrB", "TrA", "TrB", "0rD", "0rU", "1rD", "1rU", "2rD", "2rU", "3rD", "3rU", "5rD",
            "5rU", "ZrD", "ZrU",
        ],
    },
    MonoEntry {
        key: "Kdn",
        label: "Ketodeoxynononic acid (green diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::Green,
        color2: SnfgColor::Green,
        names: &["KDN"],
    },
    MonoEntry {
        key: "Neu5Ac",
        label: "N-Acetylneuraminic acid (purple diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::Purple,
        color2: SnfgColor::Purple,
        names: &[
            "SIA", "ANE5AC", "BNE5AC", "0SA", "0SB", "4SA", "4SB", "7SA", "7SB", "8SA", "8SB",
            "9SA", "9SB", "ASA", "ASB", "BSA", "BSB", "CSA", "CSB", "DSA", "DSB", "ESA", "ESB",
            "FSA", "FSB", "GSA", "GSB", "HSA", "HSB", "ISA", "ISB", "JSA", "JSB", "KSA", "KSB",
            "0sA", "0sB", "4sA", "4sB", "7sA", "7sB", "8sA", "8sB", "9sA", "9sB", "AsA", "AsB",
            "BsA", "BsB", "CsA", "CsB", "DsA", "DsB", "EsA", "EsB", "FsA", "FsB", "GsA", "GsB",
            "HsA", "HsB", "IsA", "IsB", "JsA", "JsB", "KsA",
            "KsB",
            // Not "SLB": that code collides with Gal's own GLYCAM "SLB"
            // ring-conformer variant (`Gal`'s names list, `S`+`L`+`B`),
            // and Gal is checked first.
        ],
    },
    MonoEntry {
        key: "Neu5Gc",
        label: "N-Glycolylneuraminic acid (light blue diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::LightBlue,
        color2: SnfgColor::LightBlue,
        names: &[
            "0GL", "4GL", "7GL", "8GL", "9GL", "CGL", "DGL", "EGL", "FGL", "GGL", "HGL", "IGL",
            "JGL", "KGL", "0gL", "4gL", "7gL", "8gL", "9gL", "AgL", "BgL", "CgL", "DgL", "EgL",
            "FgL", "GgL", "HgL", "IgL", "JgL", "KgL",
            "NGC", // CCD: N-glycolyl-alpha-neuraminic acid
            "NGE", // CCD: N-glycolyl-beta-neuraminic acid
        ],
    },
    MonoEntry {
        key: "Neu",
        label: "Neuraminic acid (brown diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::Brown,
        color2: SnfgColor::Brown,
        names: &["NEU"],
    },
    MonoEntry {
        key: "Sia",
        label: "Sialic acid, type unspecified (red diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::Red,
        color2: SnfgColor::Red,
        // No CCD code: the generic symbol for any of the >50 known sialic
        // acid forms, per notes.pdf footnote 6. "SIA" the CCD code is a
        // specific compound (Neu5Ac), not this generic symbol.
        names: &[],
    },
    MonoEntry {
        key: "Pse",
        label: "Pseudaminic acid (green flat diamond)",
        shape: Shape::FlatDiamond,
        color1: SnfgColor::Green,
        color2: SnfgColor::Green,
        // 6PZ's systematic name matches Pse's (Table 3) up to N-acetylation
        // at C5/C7, the biologically observed form (as SIA/Neu5Ac already
        // is Neu's own N-acetylated form).
        names: &["6PZ"],
    },
    MonoEntry {
        key: "Leg",
        label: "Legionaminic acid (yellow flat diamond)",
        shape: Shape::FlatDiamond,
        color1: SnfgColor::Yellow,
        color2: SnfgColor::Yellow,
        names: &[],
    },
    MonoEntry {
        key: "Aci",
        label: "Acinetaminic acid (pink flat diamond)",
        shape: Shape::FlatDiamond,
        color1: SnfgColor::Pink,
        color2: SnfgColor::Pink,
        names: &[],
    },
    MonoEntry {
        key: "4eLeg",
        label: "4-Epilegionaminic acid (light blue flat diamond)",
        shape: Shape::FlatDiamond,
        color1: SnfgColor::LightBlue,
        color2: SnfgColor::LightBlue,
        names: &[],
    },
    MonoEntry {
        key: "Bac",
        label: "Bacillosamine (blue hexagon)",
        shape: Shape::Hexagon,
        color1: SnfgColor::Blue,
        color2: SnfgColor::Blue,
        names: &["BAC", "0BC", "3BC", "0bC", "3bC"],
    },
    MonoEntry {
        key: "LDManHep",
        label: "L-glycero-D-manno-Heptose (green hexagon)",
        shape: Shape::Hexagon,
        color1: SnfgColor::Green,
        color2: SnfgColor::Green,
        names: &["GMH"],
    },
    MonoEntry {
        key: "Kdo",
        label: "Ketodeoxyoctonic acid (yellow hexagon)",
        shape: Shape::Hexagon,
        color1: SnfgColor::Yellow,
        color2: SnfgColor::Yellow,
        names: &["KDO"],
    },
    MonoEntry {
        key: "Dha",
        label: "3-Deoxy-lyxo-heptulosaric acid (orange hexagon)",
        shape: Shape::Hexagon,
        color1: SnfgColor::Orange,
        color2: SnfgColor::Orange,
        names: &["DHA"],
    },
    MonoEntry {
        key: "Mur",
        label: "Muramic acid (brown hexagon)",
        shape: Shape::Hexagon,
        color1: SnfgColor::Brown,
        color2: SnfgColor::Brown,
        names: &["MUR"],
    },
    MonoEntry {
        key: "DDManHep",
        label: "D-glycero-D-manno-Heptose (pink hexagon)",
        shape: Shape::Hexagon,
        color1: SnfgColor::Pink,
        color2: SnfgColor::Pink,
        names: &["289"], // CCD: D-glycero-alpha-D-manno-heptopyranose
    },
    MonoEntry {
        key: "MurNAc",
        label: "N-Acetylmuramic acid (purple hexagon)",
        shape: Shape::Hexagon,
        color1: SnfgColor::Purple,
        color2: SnfgColor::Purple,
        names: &["MUB", "AMU"], // CCD: N-acetyl-alpha/beta-muramic acid
    },
    MonoEntry {
        key: "MurNGc",
        label: "N-Glycolylmuramic acid (light blue hexagon)",
        shape: Shape::Hexagon,
        color1: SnfgColor::LightBlue,
        color2: SnfgColor::LightBlue,
        names: &[],
    },
    MonoEntry {
        key: "Api",
        label: "Apiose (blue pentagon)",
        shape: Shape::Pentagon,
        color1: SnfgColor::Blue,
        color2: SnfgColor::Blue,
        names: &["API"],
    },
    MonoEntry {
        key: "Fruc",
        label: "Fructose (green pentagon)",
        shape: Shape::Pentagon,
        color1: SnfgColor::Green,
        color2: SnfgColor::Green,
        names: &[
            "FRU", "AFRU", "BFRU", "0CA", "0CB", "1CA", "1CB", "2CA", "2CB", "3CA", "3CB", "4CA",
            "4CB", "5CA", "5CB", "WCA", "WCB", "0CD", "0CU", "1CD", "1CU", "2CD", "2CU", "3CD",
            "3CU", "4CD", "4CU", "6CD", "6CU", "WCD", "WCU", "VCD", "VCU", "UCD", "UCU", "QCD",
            "QCU", "0cA", "0cB", "1cA", "1cB", "2cA", "2cB", "3cA", "3cB", "4cA", "4cB", "5cA",
            "5cB", "WcA", "WcB", "0cD", "0cU", "1cD", "1cU", "2cD", "2cU", "3cD", "3cU", "4cD",
            "4cU", "6cD", "6cU", "WcD", "WcU", "VcD", "VcU", "UcD", "UcU", "QcD", "QcU",
        ],
    },
    MonoEntry {
        key: "Tag",
        label: "Tagatose (yellow pentagon)",
        shape: Shape::Pentagon,
        color1: SnfgColor::Yellow,
        color2: SnfgColor::Yellow,
        names: &[
            "TAG", "0JA", "0JB", "1JA", "1JB", "2JA", "2JB", "3JA", "3JB", "4JA", "4JB", "5JA",
            "5JB", "WJA", "WJB", "0JD", "0JU", "1JD", "1JU", "2JD", "2JU", "3JD", "3JU", "4JD",
            "4JU", "6JD", "6JU", "WJD", "WJU", "VJD", "VJU", "UJD", "UJU", "QJD", "QJU", "0jA",
            "0jB", "1jA", "1jB", "2jA", "2jB", "3jA", "3jB", "4jA", "4jB", "5jA", "5jB", "WjA",
            "WjB", "0jD", "0jU", "1jD", "1jU", "2jD", "2jU", "3jD", "3jU", "4jD", "4jU", "6jD",
            "6jU", "WjD", "WjU", "VjD", "VjU", "UjD", "UjU", "QjD", "QjU",
        ],
    },
    MonoEntry {
        key: "Sor",
        label: "Sorbose (orange pentagon)",
        shape: Shape::Pentagon,
        color1: SnfgColor::Orange,
        color2: SnfgColor::Orange,
        names: &[
            "SOR", "0BA", "0BB", "1BA", "1BB", "2BA", "2BB", "3BA", "3BB", "4BA", "4BB", "5BA",
            "5BB", "WBA", "WBB", "0BD", "0BU", "1BD", "1BU", "2BD", "2BU", "3BD", "3BU", "4BD",
            "4BU", "6BD", "6BU", "WBD", "WBU", "VBD", "VBU", "UBD", "UBU", "QBD", "QBU", "0bA",
            "0bB", "1bA", "1bB", "2bA", "2bB", "3bA", "3bB", "4bA", "4bB", "5bA", "5bB", "WbA",
            "WbB", "0bD", "0bU", "1bD", "1bU", "2bD", "2bU", "3bD", "3bU", "4bD", "4bU", "6bD",
            "6bU", "WbD", "WbU", "VbD", "VbU", "UbD", "UbU", "QbD", "QbU",
        ],
    },
    MonoEntry {
        key: "Psi",
        label: "Psicose (pink pentagon)",
        shape: Shape::Pentagon,
        color1: SnfgColor::Pink,
        color2: SnfgColor::Pink,
        names: &[
            "PSI", "0PA", "0PB", "1PA", "1PB", "2PA", "2PB", "3PA", "3PB", "4PA", "4PB", "5PA",
            "5PB", "WPA", "WPB", "0PD", "0PU", "1PD", "1PU", "2PD", "2PU", "3PD", "3PU", "4PD",
            "4PU", "6PD", "6PU", "WPD", "WPU", "VPD", "VPU", "UPD", "UPU", "QPD", "QPU", "0pA",
            "0pB", "1pA", "1pB", "2pA", "2pB", "3pA", "3pB", "4pA", "4pB", "5pA", "5pB", "WpA",
            "WpB", "0pD", "0pU", "1pD", "1pU", "2pD", "2pU", "3pD", "3pU", "4pD", "4pU", "6pD",
            "6pU", "WpD", "WpU", "VpD", "VpU", "UpD", "UpU", "QpD", "QpU",
        ],
    },
    // The 12 row-generic symbols (Table 1's "White" column): a residue
    // whose exact monosaccharide is unknown but whose class is not, per
    // notes.pdf footnote 10. Never resname-detected (`names` is empty);
    // present so every SNFG shape+color is representable and testable.
    MonoEntry {
        key: "Hexose",
        label: "Hexose, unspecified (white sphere)",
        shape: Shape::Sphere,
        color1: SnfgColor::White,
        color2: SnfgColor::White,
        names: &[],
    },
    MonoEntry {
        key: "HexNAc",
        label: "HexNAc, unspecified (white cube)",
        shape: Shape::Cube,
        color1: SnfgColor::White,
        color2: SnfgColor::White,
        names: &[],
    },
    MonoEntry {
        key: "Hexosamine",
        label: "Hexosamine, unspecified (white crossed cube)",
        shape: Shape::Cube,
        color1: SnfgColor::White,
        color2: SnfgColor::White,
        names: &[],
    },
    MonoEntry {
        key: "Hexuronate",
        label: "Hexuronate, unspecified (white divided diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::White,
        color2: SnfgColor::White,
        names: &[],
    },
    MonoEntry {
        key: "Deoxyhexose",
        label: "Deoxyhexose, unspecified (white cone)",
        shape: Shape::Cone,
        color1: SnfgColor::White,
        color2: SnfgColor::White,
        names: &[],
    },
    MonoEntry {
        key: "DeoxyhexNAc",
        label: "DeoxyhexNAc, unspecified (white divided cone)",
        shape: Shape::Cone,
        color1: SnfgColor::White,
        color2: SnfgColor::White,
        names: &[],
    },
    MonoEntry {
        key: "Dideoxyhexose",
        label: "Di-deoxyhexose, unspecified (white flat rectangle)",
        shape: Shape::Rectangle,
        color1: SnfgColor::White,
        color2: SnfgColor::White,
        names: &[],
    },
    MonoEntry {
        key: "Pentose",
        label: "Pentose, unspecified (white star)",
        shape: Shape::Star,
        color1: SnfgColor::White,
        color2: SnfgColor::White,
        names: &[],
    },
    MonoEntry {
        key: "Deoxynonulosonate",
        label: "Deoxynonulosonate, unspecified (white diamond)",
        shape: Shape::Diamond,
        color1: SnfgColor::White,
        color2: SnfgColor::White,
        names: &[],
    },
    MonoEntry {
        key: "Dideoxynonulosonate",
        label: "Di-deoxynonulosonate, unspecified (white flat diamond)",
        shape: Shape::FlatDiamond,
        color1: SnfgColor::White,
        color2: SnfgColor::White,
        names: &[],
    },
    MonoEntry {
        key: "Unknown",
        label: "Unknown saccharide (white flat hexagon)",
        shape: Shape::Hexagon,
        color1: SnfgColor::White,
        color2: SnfgColor::White,
        names: &[],
    },
    MonoEntry {
        key: "Assigned",
        label: "Assigned saccharide, unspecified (white pentagon)",
        shape: Shape::Pentagon,
        color1: SnfgColor::White,
        color2: SnfgColor::White,
        names: &[],
    },
];

// 75 named + 12 generic SNFG entries (87 total; Table 1/3, Neelamegham
// et al. 2019), 1587 residue-name variants (PDB CCD / common, CHARMM, GLYCAM).

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glycan::lookup;

    /// Every abbreviation in Table 3 (Neelamegham et al. 2019) plus the
    /// 12 row-generic symbols of Table 1: the full SNFG monosaccharide
    /// vocabulary this table must cover, keyed exactly as `MonoEntry::key`.
    const SPEC_ABBREVS: &[&str] = &[
        "4eLeg",
        "6dAlt",
        "6dAltNAc",
        "6dGul",
        "6dTal",
        "6dTalNAc",
        "Abe",
        "Aci",
        "All",
        "AllA",
        "AllN",
        "AllNAc",
        "Alt",
        "AltA",
        "AltN",
        "AltNAc",
        "Api",
        "Ara",
        "Bac",
        "Col",
        "DDManHep",
        "Dha",
        "Dig",
        "Fruc",
        "Fuc",
        "FucNAc",
        "Gal",
        "GalA",
        "GalN",
        "GalNAc",
        "Glc",
        "GlcA",
        "GlcN",
        "GlcNAc",
        "Gul",
        "GulA",
        "GulN",
        "GulNAc",
        "Ido",
        "IdoA",
        "IdoN",
        "IdoNAc",
        "Kdn",
        "Kdo",
        "Leg",
        "LDManHep",
        "Lyx",
        "Man",
        "ManA",
        "ManN",
        "ManNAc",
        "Mur",
        "MurNAc",
        "MurNGc",
        "Neu",
        "Neu5Ac",
        "Neu5Gc",
        "Oli",
        "Par",
        "Pse",
        "Psi",
        "Qui",
        "QuiNAc",
        "Rha",
        "RhaNAc",
        "Rib",
        "Sia",
        "Sor",
        "Tag",
        "Tal",
        "TalA",
        "TalN",
        "TalNAc",
        "Tyv",
        "Xyl",
        // Table 1's 12 white row-generics.
        "Hexose",
        "HexNAc",
        "Hexosamine",
        "Hexuronate",
        "Deoxyhexose",
        "DeoxyhexNAc",
        "Dideoxyhexose",
        "Pentose",
        "Deoxynonulosonate",
        "Dideoxynonulosonate",
        "Unknown",
        "Assigned",
    ];

    #[test]
    fn every_spec_abbreviation_has_exactly_one_table_entry() {
        for &abbrev in SPEC_ABBREVS {
            let n = MONO_TABLE.iter().filter(|e| e.key == abbrev).count();
            assert_eq!(
                n, 1,
                "{abbrev}: expected exactly one MONO_TABLE entry, found {n}"
            );
        }
        assert_eq!(
            MONO_TABLE.len(),
            SPEC_ABBREVS.len(),
            "MONO_TABLE has entries not in SPEC_ABBREVS, or vice versa"
        );
    }

    #[test]
    fn every_table_entry_has_a_nonempty_label() {
        for e in MONO_TABLE {
            assert!(!e.label.is_empty(), "{}: empty label", e.key);
        }
    }

    /// Spot checks against the spec (Table 1/2, Neelamegham et al. 2019).
    #[test]
    fn spec_spot_checks() {
        let cases: &[(&str, Shape, SnfgColor, SnfgColor)] = &[
            ("GlcNAc", Shape::Cube, SnfgColor::Blue, SnfgColor::Blue),
            ("Man", Shape::Sphere, SnfgColor::Green, SnfgColor::Green),
            ("Gal", Shape::Sphere, SnfgColor::Yellow, SnfgColor::Yellow),
            ("Fuc", Shape::Cone, SnfgColor::Red, SnfgColor::Red),
            (
                "Neu5Ac",
                Shape::Diamond,
                SnfgColor::Purple,
                SnfgColor::Purple,
            ),
            ("Xyl", Shape::Star, SnfgColor::Orange, SnfgColor::Orange),
            ("GlcA", Shape::Diamond, SnfgColor::Blue, SnfgColor::White),
            (
                "Pse",
                Shape::FlatDiamond,
                SnfgColor::Green,
                SnfgColor::Green,
            ),
            ("QuiNAc", Shape::Cone, SnfgColor::White, SnfgColor::Blue),
            // Deoxyhexose row: colors follow the hexose row's own
            // Gul/Alt/Tal columns (pdftotext column match), not a plain
            // sequential fill after Qui=Blue, Rha=Green.
            ("6dGul", Shape::Cone, SnfgColor::Orange, SnfgColor::Orange),
            (
                "6dTalNAc",
                Shape::Cone,
                SnfgColor::White,
                SnfgColor::LightBlue,
            ),
        ];
        for &(key, shape, c1, c2) in cases {
            let e = MONO_TABLE.iter().find(|e| e.key == key).unwrap();
            assert_eq!(e.shape, shape, "{key}: shape");
            assert_eq!(e.color1, c1, "{key}: color1");
            assert_eq!(e.color2, c2, "{key}: color2");
        }
    }

    /// CCD residue codes this audit added, verified individually against
    /// `data.rcsb.org/rest/v1/core/chemcomp/<CODE>`.
    #[test]
    fn slb_stays_glycam_gal_not_ccd_neu5ac() {
        // "SLB" is ambiguous: Gal's own GLYCAM ring-conformer code, and
        // separately the CCD code for beta-Neu5Ac. Gal is listed first
        // and must keep winning (see the comment on Neu5Ac's names).
        assert_eq!(lookup("SLB").unwrap().key, "Gal");
    }

    #[test]
    fn newly_mapped_ccd_codes_resolve_to_the_right_sugar() {
        let cases: &[(&str, &str)] = &[
            ("NGC", "Neu5Gc"),
            ("NGE", "Neu5Gc"),
            ("IDR", "IdoA"),
            ("BDP", "GlcA"),
            ("GTR", "GalA"),
            ("A2G", "GalNAc"),
            ("PA1", "GlcN"),
            ("RM4", "Rha"),
            ("ARB", "Ara"),
            ("MUB", "MurNAc"),
            ("AMU", "MurNAc"),
            ("X6X", "GalN"),
            ("289", "DDManHep"),
        ];
        for &(code, key) in cases {
            let e = lookup(code).unwrap_or_else(|| panic!("{code} not recognized"));
            assert_eq!(e.key, key, "{code}");
        }
    }
}
